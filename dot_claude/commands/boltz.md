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
- `--recycling N` - Recycling steps (default: 3)
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

### 5. Install Boltz (if needed)

**IMPORTANT**: Use `uv` for faster installation (not pip):

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
if ! command -v boltz &> /dev/null; then
    # Install uv first if not present
    if ! command -v uv &> /dev/null; then
        curl -LsSf https://astral.sh/uv/install.sh | sh
    fi
    source ~/.local/bin/env
    uv pip install boltz --system
fi
export PATH=\$HOME/.local/bin:\$PATH
boltz --version
"
```

**Skip this step** if using `/reserve-gpu --image boltz/boltz` (Boltz pre-installed).

### 6. Run Boltz Prediction

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:\$LD_LIBRARY_PATH
export PATH=\$HOME/.local/bin:\$PATH
cd /tmp
boltz predict input.yaml \
    --out_dir output \
    --use_msa_server \
    --recycling_steps 3 \
    --diffusion_samples 1 \
    --num_workers 0 \
    --no_kernels
"
```

**Required flags**:
- `--no_kernels`: **ESSENTIAL** - Disables cuequivariance_torch kernels (not installed by default)
- `--num_workers 0`: Avoids shared memory / bus error issues in containers
- `--use_msa_server`: Uses ColabFold MSA server (no local database needed)

**Optional flags**:
- `--recycling_steps 3`: Structure refinement iterations (default 3)
- `--diffusion_samples N`: Number of output models (use 1 for speed, 5 for production)

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
| Metric | Good Value | Description |
|--------|------------|-------------|
| `confidence_score` | > 0.8 | Overall prediction quality |
| `ptm` | > 0.8 | Predicted TM-score (protein quality) |
| `complex_plddt` | > 0.8 | Per-residue confidence |
| `ligand_iptm` | > 0.8 | Ligand positioning confidence (only for protein-ligand) |

**Example good results**:

Protein-only (ubiquitin):
```json
{"confidence_score": 0.933, "ptm": 0.918, "complex_plddt": 0.937}
```

Covalent ligands (CDK2-Cys118):
| Fragment | Warhead | Confidence | Ligand iPTM |
|----------|---------|------------|-------------|
| Frag1 | Acrylamide (Michael) | 0.92 | 0.92 |
| Frag2 | Chloroacetamide (SN2) | 0.93 | 0.94 |

## Troubleshooting

### Common Errors

**ModuleNotFoundError: cuequivariance_torch**: Add `--no_kernels` flag. This is the most common error!

**KeyError: 'name'**: Atom properties not preserved in pickle. Ensure `Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)` is called before creating molecules.

**RuntimeError: Class values must be smaller than num_classes**: Atom names contain invalid characters. Use UPPERCASE only: `CL1` not `Cl1`.

**Exit code 137 (OOM)**: Boltz downloads ~6-10GB model weights. Ensure sufficient disk space.

**DataLoader bus error**: Add `--num_workers 0` flag.

**boltz: command not found after uv install**: Add `export PATH=$HOME/.local/bin:$PATH` to the command.

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
