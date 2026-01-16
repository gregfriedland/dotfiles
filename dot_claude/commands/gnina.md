# GNINA Molecular Docking

Run GNINA deep learning molecular docking on a GPU Kubernetes node.

## Usage

```
/gnina <receptor> <ligand> [options]
```

### Arguments

**Receptor** (required): Protein structure file
- PDB file: `receptor.pdb`
- PDBQT file: `receptor.pdbqt`

**Ligand** (required): Ligand molecule(s)
- SDF file: `ligand.sdf`
- MOL2 file: `ligand.mol2`
- SMILES: `"CCO"` (will be converted)

### Options

- `--autobox <ligand.sdf>` - Reference ligand to define search box (recommended)
- `--center <x,y,z>` - Manual box center coordinates
- `--size <x,y,z>` - Box size in Angstroms (default: 25,25,25)
- `--exhaustiveness N` - Search thoroughness (default: 128)
- `--num_modes N` - Number of output poses (default: 9)
- `--covalent <chain:res:atom> <smarts>` - Covalent docking (see below)
- `--flex <chain:res,...>` - Flexible receptor residues
- `--cnn <model>` - CNN model: default, fast, dense (default: default ensemble)
- `--output <file.sdf>` - Output file name

## Examples

```bash
# Standard docking with reference ligand for box
/gnina receptor.pdb ligand.sdf --autobox crystal_ligand.sdf

# Dock SMILES to protein
/gnina receptor.pdb "CC(=O)Nc1ccccc1" --autobox ref.sdf

# High-throughput screening (fast CNN)
/gnina receptor.pdb library.sdf --autobox ref.sdf --cnn fast

# Covalent docking to Cys118 (acrylamide warhead)
/gnina receptor.pdb ligand.sdf --autobox ref.sdf --covalent A:118:SG "[CH2]=[CH]"

# Covalent docking to Cys (chloroacetamide warhead)
/gnina receptor.pdb ligand.sdf --autobox ref.sdf --covalent A:118:SG "[CH2]Cl"

# Flexible receptor residues
/gnina receptor.pdb ligand.sdf --autobox ref.sdf --flex A:118,A:119,A:120

# Manual box definition
/gnina receptor.pdb ligand.sdf --center 10.5,20.3,15.2 --size 20,20,20
```

## Instructions

When this command is invoked, perform the following steps:

### 1. Check/Setup GPU Pod

First, check if a GPU pod is already running:

```bash
kubectl get pod -n development | grep gpu-
```

If no pod exists or it's not running, invoke the `/reserve-gpu` skill:
```bash
/reserve-gpu 2
```

**Installing GNINA**: The base CUDA image requires installing CUDA runtime libs first:

```bash
# Install PyTorch to get CUDA runtime libraries
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install nvidia-cusparse-cu12 nvidia-nvtx-cu12

# Download gnina binary
# Use v1.1 for covalent docking (v1.3 has bugs with covalent mode)
wget -q https://github.com/gnina/gnina/releases/download/v1.1/gnina -O /usr/local/bin/gnina
chmod +x /usr/local/bin/gnina

# Set library path
export LD_LIBRARY_PATH=/usr/local/lib/python3.10/dist-packages/nvidia/cuda_runtime/lib:/usr/local/lib/python3.10/dist-packages/nvidia/cublas/lib:/usr/local/lib/python3.10/dist-packages/nvidia/cusparse/lib:/usr/local/lib/python3.10/dist-packages/nvidia/nvtx/lib:$LD_LIBRARY_PATH
```

### 2. Prepare Input Files

Copy receptor and ligand files to the pod:

```bash
kubectl cp receptor.pdb development/reserved-l4-gpu:/tmp/receptor.pdb
kubectl cp ligand.sdf development/reserved-l4-gpu:/tmp/ligand.sdf
kubectl cp reference.sdf development/reserved-l4-gpu:/tmp/reference.sdf  # if using autobox
```

#### Converting SMILES to SDF

If ligand is provided as SMILES, convert it first:

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
cat > /tmp/smiles_to_sdf.py << 'PYSCRIPT'
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

smiles = sys.argv[1]
output = sys.argv[2]

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

writer = Chem.SDWriter(output)
writer.write(mol)
writer.close()
print(f'Wrote {output}')
PYSCRIPT
python /tmp/smiles_to_sdf.py 'SMILES_HERE' /tmp/ligand.sdf
"
```

### 3. Run GNINA Docking

#### Standard Non-Covalent Docking

```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
cd /tmp
gnina \
    -r receptor.pdb \
    -l ligand.sdf \
    --autobox_ligand reference.sdf \
    --cnn_scoring rescore \
    --minimize \
    --exhaustiveness 128 \
    --num_modes 9 \
    -o docked.sdf
"
```

#### Covalent Docking

**CRITICAL**: GNINA expects the **POST-REACTION** (covalently bound) form of the ligand, NOT the pre-reaction form! The ligand must already have the warhead in its "reacted" state.

**Examples of SMILES transformation**:
```python
# Acrylamide (Michael addition): vinyl becomes saturated
# Pre-reaction:  C=CC(=O)N...  (has C=C)
# Post-reaction: CCC(=O)N...   (saturated, C1 bonds to Cys-SG)

# Chloroacetamide (SN2): Cl is replaced by bond to Cys
# Pre-reaction:  ...C(=O)CCl
# Post-reaction: ...C(=O)C     (Cl removed, terminal C bonds to Cys-SG)
```

```bash
kubectl exec gpu-gnina-xxx -n development -- bash -c "
cd /tmp
gnina \
    -r receptor.pdb \
    -l ligand_post_reaction.sdf \
    --autobox_ligand reference.sdf \
    --covalent_rec_atom A:118:SG \
    --covalent_lig_atom_pattern '[CH3]' \
    --covalent_optimize_lig \
    --cnn_scoring rescore \
    --minimize \
    --exhaustiveness 128 \
    --num_modes 9 \
    -o docked_covalent.sdf
"
```

**Covalent docking flags**:
| Flag | Description |
|------|-------------|
| `--covalent_rec_atom` | Protein atom: `chain:resnum:atomname` (e.g., `A:118:SG` for Cys sulfur) |
| `--covalent_lig_atom_pattern` | SMARTS for ligand attachment atom (in POST-REACTION form) |
| `--covalent_optimize_lig` | Optimize ligand-residue geometry with UFF |
| `--covalent_bond_order` | Bond order (default: 1) |

**Warhead SMARTS patterns (POST-REACTION form)**:
| Warhead | Pre-reaction | Post-reaction SMARTS | Notes |
|---------|--------------|---------------------|-------|
| Acrylamide | `C=C` | `[CH2][CH2]` or `C` | Saturated beta carbon |
| Chloroacetamide | `CCl` | `[CH2]` or `C` | Terminal C (Cl removed) |
| Epoxide | `C1OC1` | `[CH2]O` | Ring-opened |

### 4. Copy Results Back

```bash
kubectl cp development/reserved-l4-gpu:/tmp/docked.sdf ./gnina_output/docked.sdf
```

### 5. Report Results

Parse the output SDF file for scores. GNINA adds these properties to each pose:

| Property | Description | Good Value |
|----------|-------------|------------|
| `CNNscore` | CNN binding probability | > 0.7 |
| `CNNaffinity` | Predicted pKd | > 6 (nM binder) |
| `minimizedAffinity` | Vina-like score after minimization | < -7 kcal/mol |

```bash
# View scores from output
kubectl exec reserved-l4-gpu -n development -- bash -c "
grep -E '(CNNscore|CNNaffinity|minimizedAffinity)' /tmp/docked.sdf | head -20
"
```

## Default Settings

These defaults are optimized for accuracy:

| Setting | Value | Rationale |
|---------|-------|-----------|
| `--cnn_scoring rescore` | GPU CNN scoring | Best accuracy with CNN ensemble |
| `--minimize` | Enabled | Refines poses after docking |
| `--exhaustiveness 128` | High | Thorough search (16x default) |
| `--num_modes 9` | 9 poses | Standard output |

### For High-Throughput Screening

Use faster settings for large libraries:

```bash
gnina ... --cnn fast --exhaustiveness 32 --num_modes 3
```

## CNN Scoring Options

| Option | Speed | Accuracy | Use Case |
|--------|-------|----------|----------|
| `--cnn default` | Medium | Best | Standard docking |
| `--cnn fast` | Fast | Good | Virtual screening |
| `--cnn dense` | Slow | Best | Final ranking |
| `--cnn_scoring none` | Fastest | Vina only | Quick filtering |
| `--cnn_scoring refinement` | Slow | Best | Pose optimization |

## Flexible Receptor Docking

Allow specific residues to move during docking:

```bash
gnina -r receptor.pdb -l ligand.sdf \
    --autobox_ligand ref.sdf \
    --flexres A:118,A:119,A:120 \
    --minimize \
    --exhaustiveness 128 \
    -o docked_flex.sdf
```

**Note**: Flexible residues significantly increase runtime.

## Troubleshooting

### Common Errors

**No GPU detected**: Ensure `--device 0` is set or the container has GPU access.

**Out of memory**: Reduce `--exhaustiveness` or `--num_modes`.

**No poses found**:
- Check box size covers binding site
- Verify ligand is valid (try `obabel -isdf ligand.sdf -osdf`)
- Increase exhaustiveness

**Covalent docking fails**:
- **Use gnina v1.1** - v1.3 has internal errors with covalent mode
- Verify SMARTS pattern matches ligand atoms
- Check receptor atom name is correct (use `grep " SG " receptor.pdb`)
- Ensure ligand contains the reactive warhead
- Remove explicit hydrogens from ligand SDF
- Remove CONECT records from PDB: `grep -v "^CONECT" receptor.pdb > receptor_clean.pdb`

**CNN scoring warning for covalent**: CNN models aren't calibrated for covalent docking. Use for relative ranking only, not absolute scores. Vina affinity will be positive (unfavorable) due to covalent constraints - this is expected.

### GPU Pod Issues

If the pod doesn't have gnina:
```bash
kubectl exec reserved-l4-gpu -n development -- bash -c "
wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina -O /usr/local/bin/gnina
chmod +x /usr/local/bin/gnina
gnina --version
"
```

## References

- GNINA GitHub: https://github.com/gnina/gnina
- GNINA 1.3 Paper: https://link.springer.com/article/10.1186/s13321-025-00973-x
- Docker Hub: https://hub.docker.com/r/gnina/gnina
- Covalent Docking Tutorial: https://github.com/gnina/gnina#covalent-docking
