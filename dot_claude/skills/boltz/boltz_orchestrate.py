#!/usr/bin/env python3
"""Local orchestrator for parallel Boltz-2 scoring across K8s GPU pods.

Splits input compounds into batches, distributes across N pods in parallel,
and retrieves results after each batch. Handles timeouts and errors gracefully.
Supports resume: existing results in the output CSV are skipped.

Supports multi-chain proteins via --protein_config JSON file:
    [
        {"id": "A", "sequence": "MSEQ..."},
        {"id": "B", "sequence": "QSEQ..."}
    ]
The ligand chain ID is auto-assigned as the next letter after the last protein.

Usage:
    python boltz_orchestrate.py \
        --input bitopic_hits.parquet \
        --output bitopic_boltz_results.csv \
        --pods gpu-base-k7m2p4,gpu-base-vvcgxr,gpu-base-uq6gmi \
        --batch_size 10 \
        --recycling_steps 3 \
        --protein_config proteins.json
"""

import argparse
import csv
import json
import os
import queue
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

NS = "development"
CONTEXT = "west1"
MSA_REMOTE_PATH = "/workspace/target_msa.csv"
REMOTE_SCRIPT_NAME = "boltz_remote.py"

CDK2_SEQUENCE = (
    "MANFQKVEKIGEAAYGVVYKARNKLTGEVVALKKIRLAAAAIREISLLKELNHPNIVKLLDVI"
    "HTLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLL"
    "INTEGAIKLADFGLARAFGVPAVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRAL"
    "FPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQ"
    "MLHYDPNKRISAKAALAHPFFQDVT"
)

DEFAULT_PROTEIN_CONFIG = [{"id": "A", "sequence": CDK2_SEQUENCE}]


def load_protein_config(path):
    """Load protein chain config from JSON file, or return CDK2 default."""
    if path is None:
        return DEFAULT_PROTEIN_CONFIG
    with open(path) as f:
        return json.load(f)


def ligand_chain_id(protein_config):
    """Return the next chain letter after the last protein chain."""
    last = max(p["id"] for p in protein_config)
    return chr(ord(last) + 1)


def build_yaml(protein_config, smiles, msa_paths=None):
    """Build a Boltz YAML string from protein config + ligand SMILES.

    Args:
        protein_config: list of {"id": "A", "sequence": "..."} dicts
        smiles: ligand SMILES string
        msa_paths: optional dict mapping chain id -> remote msa path
    """
    yaml_smiles = smiles.replace("\\", "\\\\")
    lig_id = ligand_chain_id(protein_config)
    lines = ["version: 1", "sequences:"]
    for p in protein_config:
        lines.append("- protein:")
        lines.append(f"    id:")
        lines.append(f"    - {p['id']}")
        lines.append(f"    sequence: {p['sequence']}")
        if msa_paths and p["id"] in msa_paths:
            lines.append(f"    msa: {msa_paths[p['id']]}")
    lines.append("- ligand:")
    lines.append(f"    id:")
    lines.append(f"    - {lig_id}")
    lines.append(f'    smiles: "{yaml_smiles}"')
    lines.append("properties:")
    lines.append("- affinity:")
    lines.append(f"    binder: {lig_id}")
    return "\n".join(lines) + "\n"

# Fixed prefix columns; score columns discovered dynamically from batch CSVs
_FIXED_PREFIX = ["compound_id", "smiles", "status", "error", "batch_idx", "pod"]
_FIXED_SUFFIX = ["protein_pdb_zstd_b64", "ligand_sdf_zstd_b64"]

print_lock = threading.Lock()


def log(msg):
    with print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# kubectl helpers
# ---------------------------------------------------------------------------

def kexec(pod, cmd, timeout=120):
    """Run command on pod, return stdout."""
    full = (
        f"kubectl --context={CONTEXT} exec {pod} -n {NS} "
        f"-- bash -c {_shellquote(cmd)}"
    )
    r = subprocess.run(full, shell=True, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        detail = (r.stdout[-500:] + "\n" + r.stderr[-500:]).strip()
        raise RuntimeError(
            f"kexec on {pod} failed (rc={r.returncode}): {detail}"
        )
    return r.stdout


def kcp_to(local, pod, remote):
    cmd = f"kubectl --context={CONTEXT} cp {local} {NS}/{pod}:{remote}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f"kcp_to {pod} failed: {r.stderr[-300:]}")


def kcp_from(pod, remote, local):
    cmd = f"kubectl --context={CONTEXT} cp {NS}/{pod}:{remote} {local}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f"kcp_from {pod} failed: {r.stderr[-300:]}")


def _shellquote(s):
    """Quote a string for bash -c '...'."""
    return "'" + s.replace("'", "'\\''") + "'"


# ---------------------------------------------------------------------------
# Pod management
# ---------------------------------------------------------------------------

def cleanup_pod(pod):
    """Kill any running boltz processes on the pod."""
    log(f"[{pod}] Cleaning up existing boltz jobs...")
    try:
        kexec(pod, "pkill -9 -f 'boltz predict' 2>/dev/null; pkill -9 -f boltz_remote 2>/dev/null; true")
    except Exception:
        pass


def setup_pod(pod, script_dir, protein_config_path=None):
    """Install uv, deps, and copy remote script + protein config to pod."""
    log(f"[{pod}] Setting up...")
    cleanup_pod(pod)

    # Ensure /workspace exists
    kexec(pod, "mkdir -p /workspace")

    # Install uv if not present
    kexec(
        pod,
        "export PATH=$HOME/.local/bin:$PATH && "
        "which uv > /dev/null 2>&1 || "
        "(curl -LsSf https://astral.sh/uv/install.sh | sh)",
        timeout=120,
    )

    # Create venv if needed, then always install deps
    kexec(
        pod,
        "export PATH=$HOME/.local/bin:$PATH && "
        "test -d /opt/boltz_env || uv venv /opt/boltz_env",
        timeout=120,
    )
    kexec(
        pod,
        "export PATH=$HOME/.local/bin:$PATH && "
        "uv pip install --python /opt/boltz_env/bin/python "
        "boltz gemmi zstandard pyarrow",
        timeout=600,
    )

    # Copy remote script
    local_script = os.path.join(script_dir, REMOTE_SCRIPT_NAME)
    kcp_to(local_script, pod, f"/workspace/{REMOTE_SCRIPT_NAME}")

    # Copy protein config if provided
    if protein_config_path:
        kcp_to(protein_config_path, pod, "/workspace/protein_config.json")

    log(f"[{pod}] Setup complete")


def precompute_msa(pod, first_smiles, protein_config):
    """Run boltz on one compound with --use_msa_server to pre-compute MSA.

    Returns a dict mapping chain_id -> remote MSA path, or a single path string.
    """
    log(f"[{pod}] Pre-computing MSA ({len(protein_config)} chain(s))...")

    # Create and copy warmup YAML
    yaml_content = build_yaml(protein_config, first_smiles)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        local_yaml = f.name
    try:
        kexec(pod, "mkdir -p /workspace/msa_warmup")
        kcp_to(local_yaml, pod, "/workspace/msa_warmup/warmup.yaml")
    finally:
        os.unlink(local_yaml)

    # Run boltz predict for MSA generation
    kexec(
        pod,
        "source /opt/boltz_env/bin/activate && "
        "boltz predict /workspace/msa_warmup/warmup.yaml "
        "--out_dir /workspace/msa_warmup/output "
        "--use_msa_server "
        "--cache ~/.boltz "
        "--recycling_steps 1 "
        "--diffusion_samples 1 "
        "--num_workers 0 "
        "--no_kernels "
        "--override 2>&1",
        timeout=1800,
    )

    # Find MSA files — Boltz stores them as <name>_<chain_idx>.csv in msa/ dir.
    # For multi-chain there is one per protein chain (warmup_0.csv, warmup_1.csv, ...).
    msa_out = kexec(
        pod,
        "find /workspace/msa_warmup/output -path '*/msa/*.csv' -maxdepth 5 -type f 2>/dev/null"
    ).strip()
    if not msa_out:
        msa_out = kexec(
            pod,
            "find /workspace/msa_warmup/output -name '*.a3m' -type f 2>/dev/null"
        ).strip()
    if not msa_out:
        raise RuntimeError("MSA warmup failed: no MSA file found")

    msa_files = sorted(msa_out.strip().split("\n"))
    # Copy MSA files to canonical locations
    msa_paths = {}
    for i, p_cfg in enumerate(protein_config):
        if i < len(msa_files):
            remote_path = f"/workspace/msa_chain_{p_cfg['id']}.csv"
            kexec(pod, f"cp {msa_files[i].strip()} {remote_path}")
            msa_paths[p_cfg["id"]] = remote_path
            log(f"[{pod}] MSA chain {p_cfg['id']}: {msa_files[i].strip()} -> {remote_path}")

    return msa_paths


def distribute_msa(source_pod, target_pods, msa_paths):
    """Copy pre-computed MSA files from source pod to all other pods."""
    if not target_pods:
        return
    for chain_id, remote_path in msa_paths.items():
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            local_msa = f.name
        try:
            kcp_from(source_pod, remote_path, local_msa)
            for pod in target_pods:
                kcp_to(local_msa, pod, remote_path)
                log(f"[{pod}] MSA chain {chain_id} distributed")
        finally:
            os.unlink(local_msa)


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def process_batch(pod, batch_df, batch_num, msa_paths, recycling_steps,
                   protein_config_path=None):
    """Send one batch to a pod, run remote script, fetch results CSV."""
    remote_input = f"/workspace/batch_{batch_num}.parquet"
    remote_output = f"/workspace/batch_{batch_num}_results.csv"

    # Write batch parquet locally and copy to pod
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
        batch_df.to_parquet(f.name, index=False)
        local_input = f.name
    try:
        kcp_to(local_input, pod, remote_input)
    finally:
        os.unlink(local_input)

    # Build remote command
    cmd = (
        f"source /opt/boltz_env/bin/activate && "
        f"timeout 3500 python /workspace/{REMOTE_SCRIPT_NAME} "
        f"--input {remote_input} "
        f"--output {remote_output} "
        f"--recycling_steps {recycling_steps}"
    )
    if protein_config_path:
        cmd += " --protein_config /workspace/protein_config.json"
    if msa_paths:
        # Pass all per-chain MSA paths as comma-separated chain_id:path pairs
        msa_arg = ",".join(f"{cid}:{p}" for cid, p in msa_paths.items())
        cmd += f" --msa_files {msa_arg}"
    # Note: stderr captured by kexec for error reporting

    # Execute with 1hr local timeout
    exec_ok = True
    try:
        stdout = kexec(pod, cmd, timeout=3600)
        log(f"[{pod}] Batch {batch_num} remote output (last 500):\n{stdout[-500:]}")
    except Exception as e:
        exec_ok = False
        log(f"[{pod}] Batch {batch_num} execution error: {e}")
        # Fall through to try fetching whatever results exist

    # Fetch results CSV
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        local_output = f.name
    try:
        kcp_from(pod, remote_output, local_output)
        with open(local_output, newline="") as f:
            reader = csv.DictReader(f)
            results = list(reader)
        # Guard against empty results (e.g., stale/empty CSV)
        if not results:
            log(f"[{pod}] Batch {batch_num}: CSV fetched but empty")
            return [
                {
                    "compound_id": row["compound_id"],
                    "smiles": row["smiles"],
                    "status": "failed",
                    "error": "empty_results_csv",
                }
                for _, row in batch_df.iterrows()
            ]
        return results
    except Exception as e:
        log(f"[{pod}] Batch {batch_num}: could not fetch results: {e}")
        # Mark all compounds as failed
        return [
            {
                "compound_id": row["compound_id"],
                "smiles": row["smiles"],
                "status": "failed",
                "error": "batch_execution_failed",
            }
            for _, row in batch_df.iterrows()
        ]
    finally:
        if os.path.exists(local_output):
            os.unlink(local_output)


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

def _build_fieldnames(results):
    """Derive CSV fieldnames dynamically: fixed prefix + discovered score cols + fixed suffix."""
    seen = []
    fixed = set(_FIXED_PREFIX + _FIXED_SUFFIX)
    for r in results:
        for k in r:
            if k not in fixed and k not in seen:
                seen.append(k)
    return _FIXED_PREFIX + seen + _FIXED_SUFFIX


def save_csv(output_path, results):
    fieldnames = _build_fieldnames(results)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)


def load_csv(output_path):
    if not os.path.exists(output_path):
        return []
    with open(output_path, newline="") as f:
        return list(csv.DictReader(f))


def save_structures(output_path, results):
    """Decode zstd+b64 structure columns and save PDB/SDF files next to the CSV."""
    import base64

    try:
        import zstandard as zstd
    except ImportError:
        log("WARNING: zstandard not installed, skipping structure extraction")
        return

    out_dir = os.path.dirname(os.path.abspath(output_path))
    decomp = zstd.ZstdDecompressor()
    n_pdb = 0
    n_sdf = 0

    for r in results:
        if r.get("status") != "success":
            continue
        cid = r["compound_id"]

        pdb_b64 = r.get("protein_pdb_zstd_b64")
        if pdb_b64:
            pdb_str = decomp.decompress(base64.b64decode(pdb_b64)).decode()
            path = os.path.join(out_dir, f"{cid}_boltz_protein.pdb")
            with open(path, "w") as f:
                f.write(pdb_str)
            n_pdb += 1

        sdf_b64 = r.get("ligand_sdf_zstd_b64")
        if sdf_b64:
            sdf_str = decomp.decompress(base64.b64decode(sdf_b64)).decode()
            path = os.path.join(out_dir, f"{cid}_boltz_ligand.sdf")
            with open(path, "w") as f:
                f.write(sdf_str)
            n_sdf += 1

    log(f"Saved {n_pdb} PDB + {n_sdf} SDF files to {out_dir}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parallel Boltz-2 scoring across K8s GPU pods"
    )
    parser.add_argument("--input", required=True, help="Input parquet")
    parser.add_argument("--output", required=True, help="Output CSV (cumulative)")
    parser.add_argument(
        "--pods", required=True, help="Comma-separated pod names"
    )
    parser.add_argument("--batch_size", type=int, default=10)
    parser.add_argument("--recycling_steps", type=int, default=3)
    parser.add_argument(
        "--protein_config",
        default=None,
        help="JSON file with protein chain config: "
        '[{"id": "A", "sequence": "..."}, ...]. Defaults to CDK2.',
    )
    parser.add_argument(
        "--no_precompute_msa",
        action="store_true",
        help="Skip MSA pre-computation; each batch queries the MSA server",
    )
    args = parser.parse_args()

    pods = [p.strip() for p in args.pods.split(",")]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    protein_config = load_protein_config(args.protein_config)
    log(f"Protein config: {len(protein_config)} chain(s) "
        f"({', '.join(p['id'] for p in protein_config)}), "
        f"ligand chain {ligand_chain_id(protein_config)}")

    # Load input
    df = pd.read_parquet(args.input)
    log(f"Loaded {len(df)} compounds from {args.input}")

    # Resume: skip only successfully completed compounds
    existing = load_csv(args.output)
    success_cids = {r["compound_id"] for r in existing if r.get("status") == "success"}
    # Keep only successful results for the cumulative CSV
    existing = [r for r in existing if r.get("status") == "success"]
    todo = df[~df["compound_id"].isin(success_cids)].reset_index(drop=True)
    if success_cids:
        log(f"Resuming: {len(success_cids)} successful, {len(todo)} remaining")
    if len(todo) == 0:
        log("Nothing to do — all compounds already processed")
        return

    # Split into batches
    batches = []
    for i in range(0, len(todo), args.batch_size):
        batch_num = i // args.batch_size + 1
        batches.append((batch_num, todo.iloc[i : i + args.batch_size]))
    log(f"{len(batches)} batches of ≤{args.batch_size} across {len(pods)} pods")

    # Phase 1: Setup pods in parallel
    log("=" * 60)
    log("PHASE 1: Setting up pods")
    with ThreadPoolExecutor(max_workers=len(pods)) as ex:
        list(ex.map(
            lambda p: setup_pod(p, script_dir, args.protein_config), pods
        ))

    # Phase 2: Pre-compute MSA
    msa_paths = None
    if not args.no_precompute_msa:
        log("=" * 60)
        log("PHASE 2: Pre-computing MSA")
        first_smiles = todo.iloc[0]["smiles"]
        msa_paths = precompute_msa(pods[0], first_smiles, protein_config)
        if len(pods) > 1:
            distribute_msa(pods[0], pods[1:], msa_paths)
    else:
        log("Skipping MSA pre-computation (--no_precompute_msa)")

    # Phase 3: Process batches across pods
    log("=" * 60)
    log("PHASE 3: Processing batches")
    all_results = list(existing)
    results_lock = threading.Lock()
    batch_queue = queue.Queue()
    for b in batches:
        batch_queue.put(b)

    def worker(pod):
        while True:
            try:
                batch_num, batch_df = batch_queue.get_nowait()
            except queue.Empty:
                return
            log(f"[{pod}] Starting batch {batch_num}/{len(batches)} "
                f"({len(batch_df)} compounds)")
            try:
                batch_results = process_batch(
                    pod, batch_df, batch_num, msa_paths, args.recycling_steps,
                    protein_config_path=args.protein_config,
                )
                # Annotate results with batch/pod metadata
                for r in batch_results:
                    r["batch_idx"] = batch_num
                    r["pod"] = pod

                with results_lock:
                    all_results.extend(batch_results)
                    save_csv(args.output, all_results)

                n_ok = sum(1 for r in batch_results if r.get("status") == "success")
                log(f"[{pod}] Batch {batch_num} done: {n_ok}/{len(batch_results)} "
                    f"success (cumulative: {len(all_results)}/{len(df)})")
            except Exception as e:
                log(f"[{pod}] Batch {batch_num} FAILED: {e}")
                # Mark as failed and save
                failed = [
                    {
                        "compound_id": row["compound_id"],
                        "smiles": row["smiles"],
                        "status": "failed",
                        "error": str(e),
                        "batch_idx": batch_num,
                        "pod": pod,
                    }
                    for _, row in batch_df.iterrows()
                ]
                with results_lock:
                    all_results.extend(failed)
                    save_csv(args.output, all_results)
            finally:
                batch_queue.task_done()

    with ThreadPoolExecutor(max_workers=len(pods)) as ex:
        list(ex.map(worker, pods))

    # Extract structure files
    save_structures(args.output, all_results)

    # Summary
    log("=" * 60)
    n_ok = sum(1 for r in all_results if r.get("status") == "success")
    n_fail = len(all_results) - n_ok
    log(f"ALL DONE: {len(all_results)} total, {n_ok} success, {n_fail} failed")
    log(f"Results: {args.output}")


if __name__ == "__main__":
    main()
