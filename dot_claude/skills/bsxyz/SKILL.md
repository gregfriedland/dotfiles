---
name: bsxyz
description: Use when the user asks to split a PDB into protein and ligand, extract binding site coordinates, compute ligand center of mass, or create a .bsxyz file.
version: 1.0.0
---

# Binding Site XYZ Extraction

Split a protein-ligand complex PDB into separate protein and ligand files, then compute the ligand center of mass as binding site coordinates.

## Input

A PDB file containing a protein-ligand complex (e.g., `TARGET_PROTEIN.pdb`).

## Output Files

| File | Description |
|------|-------------|
| `TARGET_PROTEIN_prot.pdb` | Protein-only PDB (ATOM records + waters/ions, ligand HETATM removed) |
| `TARGET_PROTEIN_prot.bsxyz` | Binding site coordinates: `[x,y,z]` (ligand center of mass) |

Both output files are written to the same directory as the input PDB.

## Procedure

### 1. Split PDB into protein and ligand

Parse PDB line-by-line using PDB column format:
- **Protein lines**: All `ATOM` records, `TER`, `END`, header/remark lines, and HETATM records for common non-ligand residues
- **Ligand lines**: All other `HETATM` records (the actual small-molecule ligand)

Non-ligand HETATM residues to keep with protein:
```
HOH, WAT, NA, CL, MG, ZN, CA, K, SO4, PO4, GOL, EDO, ACE, NME
```

Print the ligand residue name(s) found for verification.

### 2. Compute ligand center of mass

Extract x, y, z coordinates from ligand HETATM lines using PDB fixed-width columns:
- x: columns 31-38
- y: columns 39-46
- z: columns 47-54

Compute the arithmetic mean of all ligand atom coordinates (unweighted center of mass).

### 3. Write output files

- Write protein PDB lines to `_prot.pdb`
- Write center of mass to `.bsxyz` as a single line: `[x.xxx,y.yyy,z.zzz]` (3 decimal places)

## Example Python Script

```python
from pathlib import Path
import numpy as np

pdb_path = Path("TARGET_PROTEIN.pdb")
out_dir = pdb_path.parent

with open(pdb_path) as f:
    lines = f.readlines()

protein_lines = []
ligand_lines = []
non_ligand = {"HOH", "WAT", "NA", "CL", "MG", "ZN", "CA", "K", "SO4", "PO4", "GOL", "EDO", "ACE", "NME"}

for line in lines:
    record = line[:6].strip()
    if record == "HETATM":
        resname = line[17:20].strip()
        if resname in non_ligand:
            protein_lines.append(line)
        else:
            ligand_lines.append(line)
    elif record == "ATOM" or record == "TER":
        protein_lines.append(line)
    elif record == "END":
        protein_lines.append(line)
        ligand_lines.append(line)
    else:
        protein_lines.append(line)

# Print ligand residue names for verification
lig_resnames = {line[17:20].strip() for line in ligand_lines if line[:6].strip() == "HETATM"}
print(f"Ligand residues found: {lig_resnames}")

# Write protein PDB
prot_path = out_dir / (pdb_path.stem + "_prot.pdb")
with open(prot_path, "w") as f:
    f.writelines(protein_lines)

# Compute ligand center of mass
coords = []
for line in ligand_lines:
    if line[:6].strip() == "HETATM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])

coords = np.array(coords)
com = coords.mean(axis=0)

# Write binding site xyz
bs_path = out_dir / (pdb_path.stem + "_prot.bsxyz")
with open(bs_path, "w") as f:
    f.write(f"[{com[0]:.3f},{com[1]:.3f},{com[2]:.3f}]\n")
```

## Naming Convention

Given input `NAME.pdb`:
- Protein: `NAME_prot.pdb`
- Binding site: `NAME_prot.bsxyz`

If the user specifies different output names, use those instead.
