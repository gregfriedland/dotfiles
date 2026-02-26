#!/usr/bin/env python3
"""Boltz-2 batch scoring — runs on a K8s GPU pod.

Processes all compounds in the input parquet through a single `boltz predict`
call on a directory of YAMLs. Model weights are cached at ~/.boltz. If an MSA
file is provided via --msa_file, it is referenced in every YAML and the MSA
server is not queried.

Supports multi-chain proteins via --protein_config JSON file:
    [{"id": "A", "sequence": "..."}, {"id": "B", "sequence": "..."}]

Usage:
    python boltz_remote.py \
        --input /workspace/batch.parquet \
        --output /workspace/batch_results.csv \
        [--msa_file /workspace/target_msa.csv] \
        [--protein_config /workspace/protein_config.json] \
        [--recycling_steps 3]
"""

import argparse
import base64
import csv
import json
import os
import shutil
import subprocess
import time

import pandas as pd
import zstandard as zstd

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


def parse_msa_files(msa_arg):
    """Parse --msa_files arg into a dict of chain_id -> path.

    Accepts either a single path (applied to first chain for backward compat)
    or comma-separated chain_id:path pairs (e.g. "A:/workspace/msa_A.csv,B:/workspace/msa_B.csv").
    """
    if msa_arg is None:
        return None
    pairs = {}
    for part in msa_arg.split(","):
        part = part.strip()
        if ":" in part:
            cid, path = part.split(":", 1)
            pairs[cid] = path
        else:
            # Single path without chain id — backward compat
            pairs["_single"] = part
    return pairs if pairs else None


def build_yaml(protein_config, smiles, msa_paths=None):
    """Build a Boltz YAML string from protein config + ligand SMILES.

    Args:
        protein_config: list of {"id": "A", "sequence": "..."} dicts
        smiles: ligand SMILES string
        msa_paths: optional dict mapping chain_id -> msa path, or
                   {"_single": path} for backward compat (first chain only)
    """
    yaml_smiles = smiles.replace("\\", "\\\\")
    lig_id = ligand_chain_id(protein_config)
    lines = ["version: 1", "sequences:"]
    for i, p in enumerate(protein_config):
        lines.append("- protein:")
        lines.append(f"    id:")
        lines.append(f"    - {p['id']}")
        lines.append(f"    sequence: {p['sequence']}")
        if msa_paths:
            if p["id"] in msa_paths:
                lines.append(f"    msa: {msa_paths[p['id']]}")
            elif "_single" in msa_paths and i == 0:
                lines.append(f"    msa: {msa_paths['_single']}")
    lines.append("- ligand:")
    lines.append(f"    id:")
    lines.append(f"    - {lig_id}")
    lines.append(f'    smiles: "{yaml_smiles}"')
    lines.append("properties:")
    lines.append("- affinity:")
    lines.append(f"    binder: {lig_id}")
    return "\n".join(lines) + "\n"

COMPRESSOR = zstd.ZstdCompressor(level=3)

# Fixed columns; score columns discovered dynamically from JSONs
_FIXED_PREFIX = ["compound_id", "smiles", "status", "error"]
_FIXED_SUFFIX = ["protein_pdb_zstd_b64", "ligand_sdf_zstd_b64"]


def scalar_items(obj):
    """Return only top-level scalar (non-dict, non-list) key-value pairs."""
    return {k: v for k, v in obj.items() if isinstance(v, (int, float, str, bool, type(None)))}


def extract_coords_from_cif(cif_path, smiles, protein_chain_ids, lig_chain_id):
    """Extract protein PDB and ligand SDF from Boltz CIF.

    Args:
        cif_path: path to the CIF file
        smiles: ligand SMILES for bond order assignment
        protein_chain_ids: list of protein chain IDs (e.g., ["A", "B"])
        lig_chain_id: ligand chain ID (e.g., "C")
    """
    import gemmi

    doc = gemmi.cif.read(cif_path)
    st = gemmi.make_structure_from_block(doc[0])

    # Protein chains → PDB
    prot_chain_set = set(protein_chain_ids)
    prot_st = gemmi.Structure()
    prot_st.cell = st.cell
    prot_st.spacegroup_hm = st.spacegroup_hm
    for model in st:
        m = gemmi.Model(model.name)
        for chain in model:
            if chain.name in prot_chain_set:
                m.add_chain(chain)
        prot_st.add_model(m)
    protein_pdb = prot_st.make_pdb_string()

    # Ligand chain → SDF via RDKit
    ligand_sdf = None
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        lig_st = gemmi.Structure()
        lig_st.cell = st.cell
        for model in st:
            m = gemmi.Model(model.name)
            for chain in model:
                if chain.name == lig_chain_id:
                    m.add_chain(chain)
            lig_st.add_model(m)

        mol = Chem.MolFromPDBBlock(lig_st.make_pdb_string(), removeHs=False)
        if mol is not None:
            template = Chem.MolFromSmiles(smiles)
            if template is not None:
                try:
                    mol = AllChem.AssignBondOrdersFromTemplate(
                        template, Chem.RemoveHs(mol)
                    )
                    mol = Chem.AddHs(mol, addCoords=True)
                except Exception:
                    pass
            from io import StringIO

            sio = StringIO()
            w = Chem.SDWriter(sio)
            w.write(mol)
            w.close()
            ligand_sdf = sio.getvalue()
    except Exception as e:
        print(f"  WARNING: ligand SDF extraction failed: {e}")

    return protein_pdb, ligand_sdf


def find_pred_base(out_dir):
    """Locate the boltz_results_* directory inside out_dir."""
    for name in sorted(os.listdir(out_dir)):
        full = os.path.join(out_dir, name)
        if os.path.isdir(full) and name.startswith("boltz_results"):
            return full
    return None


def main():
    parser = argparse.ArgumentParser(description="Boltz-2 remote batch scoring")
    parser.add_argument("--input", required=True, help="Input parquet")
    parser.add_argument("--output", required=True, help="Output CSV")
    parser.add_argument(
        "--msa_files", default=None,
        help="Pre-computed MSA files: comma-separated chain_id:path pairs "
        "(e.g. 'A:/workspace/msa_A.csv,B:/workspace/msa_B.csv') "
        "or a single path for backward compat",
    )
    parser.add_argument(
        "--protein_config",
        default=None,
        help="JSON file with protein chain config",
    )
    parser.add_argument("--recycling_steps", type=int, default=3)
    args = parser.parse_args()

    protein_config = load_protein_config(args.protein_config)
    lig_id = ligand_chain_id(protein_config)
    protein_chain_ids = [p["id"] for p in protein_config]
    print(f"Protein config: {len(protein_config)} chain(s) "
          f"({', '.join(protein_chain_ids)}), ligand chain {lig_id}")

    df = pd.read_parquet(args.input)
    print(f"Loaded {len(df)} compounds from {args.input}")

    work_dir = "/workspace/boltz_work"
    yaml_dir = os.path.join(work_dir, "yamls")
    out_dir = os.path.join(work_dir, "output")
    # Clean both yaml and output dirs to ensure a fresh start per batch
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(yaml_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # Parse MSA file paths
    msa_paths = parse_msa_files(args.msa_files)
    if msa_paths:
        print(f"MSA paths: {msa_paths}")

    # Generate YAMLs
    cid_map = {}  # safe_cid -> (compound_id, smiles)
    for _, row in df.iterrows():
        cid = row["compound_id"]
        smiles = row["smiles"]
        safe_cid = cid.replace("-", "_")
        cid_map[safe_cid] = (cid, smiles)

        yaml_content = build_yaml(
            protein_config, smiles, msa_paths=msa_paths
        )
        with open(os.path.join(yaml_dir, f"{safe_cid}.yaml"), "w") as f:
            f.write(yaml_content)

    print(f"Generated {len(cid_map)} YAMLs (msa={'provided' if msa_paths else 'server'})")

    # Run boltz predict on the directory
    cmd = (
        f"boltz predict {yaml_dir} "
        f"--out_dir {out_dir} "
        f"--cache ~/.boltz "
        f"--recycling_steps {args.recycling_steps} "
        f"--diffusion_samples 1 "
        f"--num_workers 0 "
        f"--no_kernels"
    )
    if not msa_paths:
        cmd += " --use_msa_server"
    print(f"\nRunning: {cmd}")

    start = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    elapsed = time.time() - start
    print(f"Boltz finished in {elapsed:.1f}s, rc={result.returncode}")
    if result.stdout:
        print(f"STDOUT (last 2000):\n{result.stdout[-2000:]}")
    if result.returncode != 0:
        print(f"STDERR (last 2000):\n{result.stderr[-2000:]}")

    # Collect results — score columns discovered dynamically from JSONs
    pred_base = find_pred_base(out_dir)
    results = []
    score_keys_seen = []  # ordered list of unique score column names

    for safe_cid, (cid, smiles) in cid_map.items():
        rec = {"compound_id": cid, "smiles": smiles, "status": "failed", "error": None}

        if pred_base is None:
            rec["error"] = "no_output_dir"
            results.append(rec)
            continue

        pred_dir = os.path.join(pred_base, "predictions", safe_cid)
        conf_path = os.path.join(pred_dir, f"confidence_{safe_cid}_model_0.json")
        aff_path = os.path.join(pred_dir, f"affinity_{safe_cid}.json")
        cif_path = os.path.join(pred_dir, f"{safe_cid}_model_0.cif")

        if os.path.exists(conf_path):
            with open(conf_path) as f:
                conf = json.load(f)
            for k, v in scalar_items(conf).items():
                rec[k] = v
                if k not in score_keys_seen:
                    score_keys_seen.append(k)
            rec["status"] = "success"
            print(f"  {cid}: conf={conf.get('confidence_score'):.4f}", end="")
        else:
            rec["status"] = "no_confidence_json"

        if os.path.exists(aff_path):
            with open(aff_path) as f:
                aff = json.load(f)
            for k, v in scalar_items(aff).items():
                rec[k] = v
                if k not in score_keys_seen:
                    score_keys_seen.append(k)
            print(f"  aff={aff.get('affinity_pred_value'):.3f}"
                  f"  p_bind={aff.get('affinity_probability_binary'):.3f}", end="")
        print()  # newline after per-compound output

        if os.path.exists(cif_path):
            try:
                pdb_str, sdf_str = extract_coords_from_cif(
                    cif_path, smiles, protein_chain_ids, lig_id
                )
                if pdb_str:
                    rec["protein_pdb_zstd_b64"] = base64.b64encode(
                        COMPRESSOR.compress(pdb_str.encode("utf-8"))
                    ).decode("ascii")
                if sdf_str:
                    rec["ligand_sdf_zstd_b64"] = base64.b64encode(
                        COMPRESSOR.compress(sdf_str.encode("utf-8"))
                    ).decode("ascii")
            except Exception as e:
                rec["error"] = f"coord_extraction: {e}"

        results.append(rec)

    # Build fieldnames dynamically: fixed prefix + discovered score cols + fixed suffix
    fieldnames = _FIXED_PREFIX + score_keys_seen + _FIXED_SUFFIX
    print(f"CSV columns ({len(score_keys_seen)} score cols): {score_keys_seen}")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_ok = sum(1 for r in results if r["status"] == "success")
    print(f"\nDone: {n_ok}/{len(results)} success -> {args.output}")


if __name__ == "__main__":
    main()
