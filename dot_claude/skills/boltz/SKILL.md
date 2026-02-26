---
name: boltz
description: Run Boltz-2 protein-ligand co-folding on K8s GPU pods. Supports single predictions and batch orchestration across multiple pods.
version: 1.0.0
---

# Boltz-2 Protein-Ligand Co-folding

Run Boltz-2 structure prediction on a GPU Kubernetes node to co-fold proteins with or without ligands.

## Usage

```
/boltz <protein(s)> [ligand(s)] [options]
```

### Arguments

**Proteins** (required): Specify one or more proteins by:
- HUGO gene name: `CDK2`, `EGFR`, `BRAF`
- UniProt ID: `P24941`, `Q9Y6K9`
- Comma-separated for multiple: `CDK2,EGFR`

**Ligands** (optional): Specify one or more ligands by:
- SMILES string: `"CC(=O)Nc1ccccc1"`
- CCD code: `ATP`, `NAD`
- Comma-separated for multiple: `"CCO,ATP"`

### Options

- `--covalent <chain>:<residue>:<atom>-<ligand_chain>:<atom>` - Add covalent bond constraint
  - Example: `--covalent A:118:SG-B:C1` (Cys118 SG to ligand atom C1)
- `--samples N` - Number of diffusion samples (default: 1 for speed, increase for production)
- `--recycling N` - Recycling steps (default: 3, use 10 for production)
- `--seed N` - Random seed for reproducibility
- `--output DIR` - Output directory (default: ./boltz_output)

## Examples

```bash
# Protein-only structure prediction (no ligand)
/boltz ubiquitin
/boltz CDK2

# Simple protein-ligand docking (non-covalent)
/boltz CDK2 ATP

# Protein with SMILES ligand
/boltz EGFR "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"

# Covalent ligand (Michael addition to Cys)
/boltz CDK2 "CCC(=O)N1CCC(CC1)C(=O)NC2=CC=C(N)C(=C2)[N+]([O-])=O" --covalent A:118:SG-B:C1

# Multiple proteins (complex)
/boltz CDK2,CCNA2 ATP

# Multiple ligands
/boltz BRAF "CC(=O)Nc1ccccc1,ATP"

# Production quality (recycling=10)
/boltz CDK2 ATP --recycling 10

# Multiple diffusion samples for conformational diversity
/boltz CDK2 ATP --samples 10

# Ternary protein complex with ligand
/boltz SKP1,SKP2,CKS1B "O=C1CCCN1CCCNC1=NC=NC2=C1SC(C1=CC=C(F)C=C1)=C2" --recycling 10

# Reproducible run with fixed seed
/boltz CDK2 ATP --seed 42
```

## Instructions

When this command is invoked, perform the following steps:

### 1. Check/Setup GPU Pod

First, check if a GPU pod is already running:

```bash
kubectl get pod reserved-l4-gpu -n development
```

If no pod exists or it's not running, invoke the `/reserve-gpu` skill to create one:
- Default to 2 hours for Boltz runs (model download + prediction)
- **Recommended**: Use `/reserve-gpu 2 --image boltz/boltz` for pre-installed Boltz
- Wait for the pod to be in Running state

### 2. Resolve Protein Sequences

For each protein specified:

**If HUGO gene name** (alphabetic, no numbers at start):
```bash
kubectl exec reserved-l4-gpu -n development -- curl -s 'https://rest.uniprot.org/uniprotkb/search?query=gene_exact:GENE_NAME+AND+organism_id:9606+AND+reviewed:true&format=fasta&size=1'
```
Replace `GENE_NAME` with the actual gene name (e.g., CDK2, EGFR).

**Important**: Some genes like UBB/UBC return polyubiquitin (multiple repeats). For ubiquitin, use the canonical 76aa monomer:
```
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
```

**If UniProt ID** (starts with letter, has numbers):
```bash
kubectl exec reserved-l4-gpu -n development -- curl -s 'https://rest.uniprot.org/uniprotkb/UNIPROT_ID.fasta'
```

Store the sequences for YAML generation.

### 3. Prepare Ligand Entries (if ligands specified)

**Important**: The `--custom_ccd` CLI flag does not exist in Boltz. Custom ligands must be injected directly into `~/.boltz/mols/`.

#### For Covalent Ligands

**CRITICAL**: The SMILES must represent the ligand **AFTER** the covalent reaction (the product bonded to protein):

```python
# Michael addition: acrylamide (C=C) becomes saturated (C-C) after Cys attack
# Original: C=CC(=O)N... -> Product: CCC(=O)N...
MICHAEL_PRODUCT = "CCC(=O)N1CCC(CC1)C(=O)NC2=CC=C(N)C(=C2)[N+]([O-])=O"

# SN2: chloroacetamide loses Cl, terminal C bonds to Cys-SG
# Original: ...C(=O)CCl -> Product: ...C(=O)C
SN2_PRODUCT = "[O-][N+](=O)C1=CC(Cl)=CC=C1N2CCN(CC2)C(=O)C"
```

#### Create Custom CCD Entry Script

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
cat > /tmp/create_ccd.py << 'PYSCRIPT'
import pickle
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# CRITICAL: Enable property pickling (properties don't survive pickle by default!)
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

def create_ccd_mol(smiles, ccd_code):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles}')
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass

    # Assign atom names - MUST use UPPERCASE only!
    # Wrong: Cl1 -> causes 'Class values must be smaller than num_classes'
    # Right: CL1
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol().upper()  # Force uppercase
        if symbol not in atom_counts:
            atom_counts[symbol] = 0
        atom_counts[symbol] += 1
        atom_name = f'{symbol}{atom_counts[symbol]}'

        # Required atom properties
        atom.SetProp('name', atom_name)
        atom.SetProp('alt_name', atom_name)
        atom.SetBoolProp('leaving_atom', False)

    # Required mol properties
    mol.SetProp('MOL_NAME', ccd_code)
    mol.SetProp('_Name', ccd_code)
    return mol

def save_ccd(mol, ccd_code):
    mol_dir = Path.home() / '.boltz' / 'mols'
    mol_dir.mkdir(parents=True, exist_ok=True)
    pkl_path = mol_dir / f'{ccd_code}.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(mol, f, protocol=pickle.HIGHEST_PROTOCOL)

    # VERIFY pickle saved correctly
    with open(pkl_path, 'rb') as f:
        mol_check = pickle.load(f)
    try:
        _ = mol_check.GetAtomWithIdx(0).GetProp('name')
        print(f'Created {pkl_path} - verified OK')
    except KeyError:
        raise RuntimeError(f'Atom properties not preserved in {pkl_path}!')

    # Print heavy atom names to identify attachment point
    print('Heavy atom names (for covalent bond constraint):')
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            print(f'  idx {atom.GetIdx()}: {atom.GetSymbol()} -> {atom.GetProp(\"name\")}')

if __name__ == '__main__':
    smiles = sys.argv[1]
    ccd_code = sys.argv[2]
    mol = create_ccd_mol(smiles, ccd_code)
    save_ccd(mol, ccd_code)
PYSCRIPT
python /tmp/create_ccd.py 'SMILES_HERE' 'CCD_CODE'
"
```

Use CCD codes like `LG1`, `LG2`, `XL1` etc. for custom SMILES ligands.

#### Identify Attachment Atom for Covalent Bond

After creating the CCD entry, identify the atom that bonds to Cys-SG:

| Warhead Type | Attachment Atom | Notes |
|--------------|-----------------|-------|
| Michael addition (acrylamide) | Beta carbon (C1) | First carbon in saturated product |
| SN2 (chloroacetamide) | Terminal methyl carbon | Carbon adjacent to carbonyl |
| Epoxide | Epoxide carbon | Carbon attacked by nucleophile |

Look at the printed atom names and identify the correct one for your `--covalent` constraint.

#### For Non-Covalent SMILES Ligands

Same process, but no constraint needed in YAML.

#### For CCD Ligands (ATP, NAD, etc.)

Use the code directly - Boltz has built-in CCD entries. No custom entry needed.

### 4. Create YAML Input File

Generate the Boltz YAML configuration:

**Protein-only (no ligand):**
```yaml
version: 1
sequences:
- protein:
    id:
    - A
    sequence: PROTEIN_SEQUENCE_HERE
```

**Protein + ligand:**
```yaml
version: 1
sequences:
- protein:
    id:
    - A
    sequence: PROTEIN_SEQUENCE_HERE
- ligand:
    id:
    - B
    ccd: LIGAND_CCD_CODE  # or smiles: SMILES_STRING for simple cases
```

**With covalent constraint:**
```yaml
version: 1
sequences:
- protein:
    id:
    - A
    sequence: PROTEIN_SEQUENCE_HERE
- ligand:
    id:
    - B
    ccd: LG1
constraints:
- bond:
    atom1:
    - A        # Protein chain
    - 118      # Residue number (1-indexed)
    - SG       # Atom name (SG for Cys sulfur)
    atom2:
    - B        # Ligand chain
    - 1        # Always 1 for ligand
    - C1       # Attachment atom from CCD entry
```

For multiple proteins, assign chain IDs A, B, C, etc.
For multiple ligands, use subsequent chain IDs.

**Multi-chain protein complex (no ligand):**
```yaml
version: 1
sequences:
- protein:
    id:
    - A
    sequence: PROTEIN1_SEQUENCE
- protein:
    id:
    - B
    sequence: PROTEIN2_SEQUENCE
- protein:
    id:
    - C
    sequence: PROTEIN3_SEQUENCE
```

**Multi-chain complex + ligand:**
```yaml
version: 1
sequences:
- protein:
    id:
    - A
    sequence: PROTEIN1_SEQUENCE
- protein:
    id:
    - B
    sequence: PROTEIN2_SEQUENCE
- protein:
    id:
    - C
    sequence: PROTEIN3_SEQUENCE
- ligand:
    id:
    - D
    smiles: "SMILES_STRING"  # or ccd: CODE
```

**With interface constraints from a PDB structure:**

Interface constraints force specific protein chains to be near each other's interface residues. Extract interface residues from a reference PDB (residues within 5 A of another chain), map PDB residue numbers to Boltz 1-based UniProt sequence positions, then use `pocket` constraints:

```yaml
version: 1
sequences:
- protein:
    id: [A]
    sequence: PROTEIN1_SEQUENCE
- protein:
    id: [B]
    sequence: PROTEIN2_SEQUENCE
- protein:
    id: [C]
    sequence: PROTEIN3_SEQUENCE
constraints:
  # Chain B should be near chain A's interface residues
  - pocket:
      binder: B
      contacts: [[A, 67], [A, 81], [A, 82], ...]
      max_distance: 6
      force: true
  # Chain A should be near chain B's interface residues
  - pocket:
      binder: A
      contacts: [[B, 95], [B, 96], [B, 97], ...]
      max_distance: 6
      force: true
  # Chain C should be near chain B's interface residues
  - pocket:
      binder: C
      contacts: [[B, 164], [B, 167], ...]
      max_distance: 6
      force: true
```

**Key findings from testing interface constraints:**
- Protein-protein interface constraints provide only **marginal improvement** -- Boltz typically finds the correct complex arrangement without hints (0.72->0.59 A RMSD for SKP1-SKP2-CKS1B).
- **Ligand pocket constraints are much more impactful** -- CDK2+7TW iPTM jumped 0.873->0.979 with pocket forcing.
- For multi-chain complexes, Boltz achieves sub-angstrom RMSD to crystal structures even without constraints.

### 5. Install Boltz (if needed)

**IMPORTANT**: The `boltz/boltz` Docker image does NOT have python3, pip, or curl pre-installed. You must install system packages first, then use `uv` with a venv (system pip is externally managed on Debian 12).

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
# Step 1: Install system packages (python3, pip, venv, curl)
apt-get update -qq && apt-get install -y -qq python3 python3-pip python3-venv curl

# Step 2: Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=\$HOME/.local/bin:\$PATH

# Step 3: Create venv and install boltz (system pip is externally managed)
python3 -m venv /opt/boltz_env
. /opt/boltz_env/bin/activate
uv pip install boltz

echo BOLTZ_INSTALLED
"
```

**When running boltz after venv install**, always activate the venv first:
```bash
export PATH=$HOME/.local/bin:$PATH
. /opt/boltz_env/bin/activate
boltz predict ...
```

**Note**: `boltz --version` does not exist (exit code 2). To verify installation, just run `boltz predict --help`.

### 6. Run Boltz Prediction

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:\$LD_LIBRARY_PATH
export PATH=\$HOME/.local/bin:\$PATH
cd /tmp
boltz predict input.yaml \
    --out_dir output \
    --use_msa_server \
    --recycling_steps 10 \
    --diffusion_samples 1 \
    --num_workers 0 \
    --no_kernels
"
```

### `boltz predict` CLI Flags Reference

#### Required flags (always use these on K8s L4 pods)

| Flag | Value | Why |
|------|-------|-----|
| `--no_kernels` | (no value) | **ESSENTIAL** -- Disables cuequivariance_torch custom CUDA kernels (not installed). Without this flag you get `ModuleNotFoundError: cuequivariance_torch`. |
| `--num_workers 0` | `0` | Avoids shared memory / bus error in containers. DataLoader workers crash with SIGBUS in K8s pods due to limited `/dev/shm`. |
| `--use_msa_server` | (no value) | Uses the ColabFold MSA server (`https://api.colabfold.com`) instead of requiring a local sequence database. Without this, you need ~2TB of sequence DBs. |

#### Prediction quality flags

| Flag | Default | Recommended | Description |
|------|---------|-------------|-------------|
| `--recycling_steps N` | `3` | `10` for production | Number of structure refinement iterations. Each recycle feeds the predicted structure back through the model to refine it. Higher values improve accuracy at the cost of runtime. |
| `--diffusion_samples N` | `1` | `1` for screening, `5-10` for production | Number of independent structure samples to generate. Each sample is a separate diffusion trajectory, producing a separate model. Output files: `model_0.cif`, `model_1.cif`, ..., `model_{N-1}.cif`. |
| `--seed N` | random | set for reproducibility | Random seed for the diffusion process. Use a fixed seed (e.g., `--seed 42`) when you want reproducible results across runs. Different seeds produce different diffusion trajectories. |

#### Practical impact of recycling_steps and diffusion_samples

**Recycling steps** (structure refinement within one model):
- `--recycling_steps 3`: Fast (~30s for a 300aa protein-ligand). Good for quick screening.
- `--recycling_steps 10`: ~2-3x slower but consistently better metrics. CDK2+RZ-088049 ligand iPTM improved 0.940->0.965. CKS1B confidence improved 0.862->0.893. **Use 10 for production runs.**
- Higher values (>10) give diminishing returns.

**Diffusion samples** (independent structure predictions):
- `--diffusion_samples 1`: One model. Fast. Good for initial screening.
- `--diffusion_samples 10`: 10 models, ~6x slower than 1 sample. Useful when you want to explore conformational diversity or find the best pose. Models are ranked by confidence but may find different binding pockets.
- Each sample runs the full diffusion + recycling pipeline independently -- they are NOT refinements of each other.
- **Key finding**: 10 samples of CDK2+7TW all found the same (wrong) pocket -- increasing samples does NOT guarantee finding the crystallographic pocket. Use pocket constraints instead.

**Timing on NVIDIA L4 (24GB VRAM):**

| Complex | recycling=3, samples=1 | recycling=10, samples=1 | recycling=3, samples=10 |
|---------|----------------------|------------------------|------------------------|
| 1 protein (~300aa) + ligand | ~30-35s | ~1-2min | ~2min |
| 3 proteins (~666aa total) | ~10min | ~16min | ~60min (est) |

#### Output flags

| Flag | Default | Description |
|------|---------|-------------|
| `--out_dir DIR` | `./` | Output directory. Creates `boltz_results_<input_name>/` inside. |

#### Output structure

```
out_dir/
  boltz_results_<name>/
    predictions/<name>/
      <name>_model_0.cif          # Best predicted structure
      <name>_model_1.cif          # (if diffusion_samples > 1)
      confidence_<name>_model_0.json  # Confidence metrics
      confidence_<name>_model_1.json  # (if diffusion_samples > 1)
    msa/                             # Multiple sequence alignments
    processed/                       # Preprocessed input data
```

### 7. Copy Results Back

```bash
kubectl cp development/reserved-l4-gpu:/tmp/output ./boltz_output
```

### 8. Report Results

After completion, report:
- Location of output files
- Confidence scores from `confidence_*.json`
- Best model file path (`*_model_0.cif`)

**Confidence metrics interpretation**:

| Metric | Good | Decent | Poor | Description |
|--------|------|--------|------|-------------|
| `confidence_score` | > 0.9 | 0.8-0.9 | < 0.8 | Composite of pTM + iPTM. Overall prediction quality. |
| `ptm` | > 0.95 | 0.85-0.95 | < 0.85 | Predicted TM-score for protein fold quality. |
| `ligand_iptm` | > 0.9 | 0.8-0.9 | < 0.8 | **Most important for binding.** How confidently the ligand is placed. |
| `protein_iptm` | > 0.9 | 0.8-0.9 | < 0.8 | Protein-protein interface confidence (multi-chain only). |
| `complex_plddt` | > 0.9 | 0.8-0.9 | < 0.8 | Per-residue confidence averaged across the complex. |
| `complex_iplddt` | > 0.9 | 0.7-0.9 | < 0.7 | Per-residue confidence at the interface only. More sensitive to binding quality. |
| `complex_ipde` | < 0.8 | 0.8-1.5 | > 1.5 | Interface predicted distance error (A). **Lower is better.** |

**Key metrics relationships:**
- `ligand_iptm` -- primary metric for ligand binding confidence (pose quality, not affinity)
- `complex_iplddt` -- more sensitive than `complex_plddt` for interface quality
- `complex_ipde` -- most sensitive to wrong pocket or wrong orientation (**lower = better**)
- `confidence_score` = weighted combination of `ptm` and `iptm`

**Multi-chain `pair_chains_iptm` matrix:**
For multi-chain complexes, the confidence JSON includes a pairwise iPTM matrix showing how confidently each pair of chains is positioned relative to each other. Diagonal = intra-chain pTM. Off-diagonal = inter-chain iPTM for that pair.

```json
"pair_chains_iptm": {
    "0": {"0": 0.93, "1": 0.92, "2": 0.70},
    "1": {"0": 0.72, "1": 0.85, "2": 0.80},
    "2": {"0": 0.74, "1": 0.92, "2": 0.97}
}
```

**Example results from production runs:**

| System | Confidence | Ligand iPTM | Complex iPDE | Notes |
|--------|-----------|-------------|-------------|-------|
| CDK2 + RZ-088049 (r=10) | 0.959 | 0.965 | 0.58 | Excellent binding |
| CDK2 + 7TW (pocket forced) | 0.944 | 0.979 | 0.60 | Best with pocket constraint |
| CDK2 + 7TW (unconstrained) | 0.895 | 0.856 | 1.85 | Wrong pocket without constraint |
| SKP1-SKP2-CKS1B ternary | 0.819 | -- | 4.02 | Protein-only, 0.72 A RMSD to PDB |
| SKP1-SKP2-CKS1B + 7TW | 0.819 | 0.960 | 2.06 | 7TW at SKP2-CKS1B interface |

---

## Batch Orchestration

For scoring many compounds in parallel across multiple GPU pods, use the orchestration scripts bundled with this skill.

### Overview

The orchestration system consists of two scripts:
- **`boltz_orchestrate.py`** (runs locally): Splits compounds into batches, distributes across N pods, collects results
- **`boltz_remote.py`** (runs on each pod): Processes a batch of compounds through `boltz predict`

### Prerequisites

1. Reserve multiple GPU pods (e.g., 3 pods for parallel processing):
   ```bash
   # Reserve pods with unique names
   /reserve-k8s-pod gpu-base-1 --image boltz/boltz
   /reserve-k8s-pod gpu-base-2 --image boltz/boltz
   /reserve-k8s-pod gpu-base-3 --image boltz/boltz
   ```

2. Prepare an input parquet file with columns:
   - `compound_id` (str): Unique compound identifier
   - `smiles` (str): SMILES string for each compound

### Usage

```bash
python ~/.claude/skills/boltz/boltz_orchestrate.py \
    --input compounds.parquet \
    --output results.csv \
    --pods gpu-base-1,gpu-base-2,gpu-base-3 \
    --batch_size 10 \
    --recycling_steps 3 \
    --protein_config proteins.json
```

### Orchestration Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input` | Yes | -- | Input parquet file with `compound_id` and `smiles` columns |
| `--output` | Yes | -- | Output CSV path (cumulative, supports resume) |
| `--pods` | Yes | -- | Comma-separated list of K8s pod names |
| `--batch_size` | No | 10 | Number of compounds per batch sent to each pod |
| `--recycling_steps` | No | 3 | Boltz recycling steps (use 10 for production) |
| `--protein_config` | No | CDK2 | JSON file with protein chain config (see below) |
| `--no_precompute_msa` | No | False | Skip MSA pre-computation; each batch queries the MSA server |

### Multi-Chain Protein Config

To use a multi-chain protein target, create a JSON file:

```json
[
    {"id": "A", "sequence": "MSEQ..."},
    {"id": "B", "sequence": "QSEQ..."}
]
```

The ligand chain ID is auto-assigned as the next letter after the last protein chain (e.g., chain C for a 2-chain protein). If `--protein_config` is not provided, the default CDK2 single-chain config is used.

### How It Works

1. **Setup phase**: Installs uv, creates venv, installs `boltz gemmi zstandard pyarrow` on each pod in parallel
2. **MSA phase**: Pre-computes MSA on the first pod using one compound, then distributes the cached MSA file to all pods (avoids redundant MSA server queries)
3. **Processing phase**: Distributes compound batches across pods using a thread pool with a work-stealing queue. Each pod processes one batch at a time.
4. **Resume support**: Existing successful results in the output CSV are skipped on restart. Failed compounds are retried.

### Output Format

The output CSV contains:
- `compound_id`, `smiles`: Input compound info
- `status`: `success` or `failed`
- `error`: Error message if failed
- `batch_idx`, `pod`: Processing metadata
- Boltz confidence scores (dynamically discovered): `confidence_score`, `ptm`, `ligand_iptm`, `complex_plddt`, `complex_iplddt`, `complex_ipde`, etc.
- Boltz affinity scores (if available): `affinity_pred_value`, `affinity_probability_binary`
- `protein_pdb_zstd_b64`: Zstd-compressed, base64-encoded protein PDB
- `ligand_sdf_zstd_b64`: Zstd-compressed, base64-encoded ligand SDF

### Customization

The orchestration scripts default to CDK2 as the protein target. To change the target:

- Use `--protein_config proteins.json` with a JSON array of `{"id": "X", "sequence": "..."}` entries.
- Supports any number of protein chains (single-chain, dimers, trimers, etc.)
- The ligand chain is automatically assigned the next letter after the last protein chain.

---

## Troubleshooting

### Common Errors

**ModuleNotFoundError: cuequivariance_torch**: Add `--no_kernels` flag. This is the most common error!

**KeyError: 'name'**: Atom properties not preserved in pickle. Ensure `Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)` is called before creating molecules.

**RuntimeError: Class values must be smaller than num_classes**: Atom names contain invalid characters. Use UPPERCASE only: `CL1` not `Cl1`.

**Exit code 137 (OOM)**: Boltz downloads ~6-10GB model weights. Ensure sufficient disk space.

**DataLoader bus error**: Add `--num_workers 0` flag.

**boltz: command not found after uv install**: Add `export PATH=$HOME/.local/bin:$PATH` to the command.

**boltz/boltz image has no python3/curl/pip**: The Docker image is minimal. Run `apt-get update && apt-get install -y python3 python3-pip python3-venv curl` first.

**Externally managed Python (Debian 12)**: `uv pip install boltz --system` fails. Create a venv: `python3 -m venv /opt/boltz_env && . /opt/boltz_env/bin/activate && uv pip install boltz`.

**boltz --version exits with code 2**: This flag doesn't exist. Boltz is installed correctly -- use `boltz predict --help` to verify.

**kubectl exec i/o timeout on long runs**: The kubectl connection may drop on predictions >15 min, but the process continues running on the pod. Check results with `kubectl exec ... -- ls /tmp/output/` after the timeout.

**Pod timer expired mid-run**: If the pod's sleep timer expires, the pod terminates. Create a new pod (default 5h), reinstall boltz, and rerun. For large jobs (multiple ternary complexes), ensure enough time.

**Wrong ligand pocket found**: Boltz may consistently place a ligand in a non-crystallographic pocket. Increasing `--diffusion_samples` does NOT fix this -- all samples may find the same wrong pocket. Use `pocket` constraints with `force: true` instead.

### GPU Pod Issues

If the pod is not responding or crashed:
```bash
kubectl delete pod reserved-l4-gpu -n development
# Then re-run /reserve-gpu (optionally with --image boltz/boltz)
```

## Recommended Setup for Repeated Use

For fastest iteration, use the pre-built Boltz Docker image:
```bash
/reserve-gpu 4 --image boltz/boltz
```

This skips installation and has model weights pre-cached.

## References

- Boltz-2 GitHub: https://github.com/jwohlwend/boltz
- Boltz Docker Hub: https://hub.docker.com/r/boltz/boltz
- UniProt REST API: https://www.uniprot.org/help/api
