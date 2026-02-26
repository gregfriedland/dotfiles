# PyMOL Session Creator

Create a PyMOL session with standardized visualization for protein-ligand complexes.

## Usage

```
/pymol <pdb_files> [ligand_files] [options]
```

### Arguments

**PDB files** (required): One or more protein structure files
- Single file: `protein.pdb`
- Multiple files: `protein1.pdb protein2.pdb`
- Glob pattern: `*.pdb` or `ensemble/*.pdb`

**Ligand files** (optional): Ligand molecules to load
- SDF file: `ligand.sdf`
- Multiple poses: `docked_poses.sdf`

### Options

- `--output <file.pse>` - Output session file (default: session.pse)
- `--interface <distance>` - Interface residue distance cutoff (default: 5.0 Å)
- `--ensemble` - Treat multiple PDBs as ensemble states in one object
- `--no-align` - Skip backbone alignment
- `--color-by <scheme>` - Color scheme: chain, spectrum, ss (default: chain)

## Examples

```bash
# Single protein with docked ligand
/pymol receptor.pdb ligand.sdf

# Multiple proteins aligned
/pymol protein1.pdb protein2.pdb protein3.pdb

# Ensemble as states
/pymol backrub_*.pdb --ensemble --output ensemble.pse

# Complex with custom interface distance
/pymol complex.pdb --interface 4.0

# Protein-ligand complex from docking
/pymol receptor.pdb docked_poses.sdf --output docking_results.pse
```

## Instructions

When this command is invoked, perform the following steps:

### 1. Parse Arguments

Extract:
- `PDB_FILES`: List of PDB files (can be glob pattern)
- `LIGAND_FILES`: List of ligand files (SDF, MOL2)
- `OUTPUT`: Output .pse file path (default: session.pse)
- `INTERFACE_DIST`: Distance cutoff for interface residues (default: 5.0)
- `ENSEMBLE_MODE`: Whether to combine PDBs as states
- `DO_ALIGN`: Whether to align structures (default: True)
- `COLOR_SCHEME`: How to color proteins

### 2. Expand File Patterns

```python
import glob
pdb_files = []
for pattern in PDB_FILES:
    pdb_files.extend(sorted(glob.glob(pattern)))
```

### 3. Generate PyMOL Script

Create a Python script that PyMOL will execute:

```python
#!/usr/bin/env python3
"""PyMOL session generator script."""

import pymol
from pymol import cmd
import os
import glob

# Configuration
PDB_FILES = {PDB_FILES_LIST}  # List of PDB paths
LIGAND_FILES = {LIGAND_FILES_LIST}  # List of ligand paths
OUTPUT_FILE = "{OUTPUT_PSE}"
INTERFACE_DIST = {INTERFACE_DIST}
ENSEMBLE_MODE = {ENSEMBLE_MODE}
DO_ALIGN = {DO_ALIGN}
COLOR_SCHEME = "{COLOR_SCHEME}"

def setup_visualization():
    """Configure PyMOL settings for publication-quality images."""
    # Background and rendering
    cmd.bg_color("black")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)

    # Cartoon settings
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("cartoon_transparency", 0.0)

    # Stick settings
    cmd.set("stick_radius", 0.15)
    cmd.set("stick_ball", 0)

    # Line settings
    cmd.set("line_width", 2)

    # Label settings
    cmd.set("label_size", 20)
    cmd.set("label_font_id", 7)

def load_proteins():
    """Load protein PDB files."""
    protein_objects = []

    if ENSEMBLE_MODE and len(PDB_FILES) > 1:
        # Load as states in single object
        ensemble_name = "ensemble"
        for i, pdb_file in enumerate(PDB_FILES):
            if i == 0:
                cmd.load(pdb_file, ensemble_name)
            else:
                cmd.load(pdb_file, ensemble_name, state=i+1)
        protein_objects.append(ensemble_name)
        print(f"Loaded {len(PDB_FILES)} structures as states in '{ensemble_name}'")
    else:
        # Load as separate objects
        for pdb_file in PDB_FILES:
            obj_name = os.path.splitext(os.path.basename(pdb_file))[0]
            # Sanitize name for PyMOL
            obj_name = obj_name.replace("-", "_").replace(" ", "_")
            cmd.load(pdb_file, obj_name)
            protein_objects.append(obj_name)
            print(f"Loaded {pdb_file} as '{obj_name}'")

    return protein_objects

def load_ligands():
    """Load ligand files (SDF, MOL2)."""
    ligand_objects = []

    for lig_file in LIGAND_FILES:
        obj_name = os.path.splitext(os.path.basename(lig_file))[0]
        obj_name = "lig_" + obj_name.replace("-", "_").replace(" ", "_")

        # Load SDF - each molecule becomes a state
        cmd.load(lig_file, obj_name)
        ligand_objects.append(obj_name)

        n_states = cmd.count_states(obj_name)
        print(f"Loaded {lig_file} as '{obj_name}' ({n_states} states)")

    return ligand_objects

def align_structures(protein_objects):
    """Align all protein backbones to the first structure."""
    if len(protein_objects) < 2 or not DO_ALIGN:
        return

    reference = protein_objects[0]
    for target in protein_objects[1:]:
        # Align on CA atoms (backbone)
        result = cmd.align(f"{target} and name CA", f"{reference} and name CA")
        rmsd = result[0]
        n_atoms = result[1]
        print(f"Aligned {target} to {reference}: RMSD={rmsd:.2f} A ({n_atoms} atoms)")

def style_proteins(protein_objects):
    """Apply cartoon style to proteins."""
    for obj in protein_objects:
        # Show as cartoon
        cmd.show("cartoon", obj)
        cmd.hide("lines", obj)
        cmd.hide("sticks", obj)

        # Color by scheme
        if COLOR_SCHEME == "chain":
            cmd.util.cbc(obj)  # Color by chain
        elif COLOR_SCHEME == "spectrum":
            cmd.spectrum("count", "rainbow", obj)
        elif COLOR_SCHEME == "ss":
            cmd.color("red", f"{obj} and ss h")  # Helices
            cmd.color("yellow", f"{obj} and ss s")  # Sheets
            cmd.color("green", f"{obj} and ss l+''")  # Loops

        # Remove hydrogens from protein display
        cmd.hide("everything", f"{obj} and hydro")

def style_ligands(ligand_objects):
    """Apply stick style to ligands, remove hydrogens."""
    for obj in ligand_objects:
        # Show as sticks
        cmd.show("sticks", obj)
        cmd.hide("cartoon", obj)
        cmd.hide("lines", obj)

        # Remove hydrogens
        cmd.hide("everything", f"{obj} and hydro")

        # Color by element with carbon in unique color
        cmd.util.cbag(obj)  # Color by atom, green carbons

        # Make ligand more visible
        cmd.set("stick_radius", 0.2, obj)

def show_interface_residues(protein_objects, ligand_objects):
    """Show protein residues within INTERFACE_DIST of ligands as lines."""
    if not ligand_objects:
        return

    # Create selection for all ligands
    all_ligands = " or ".join(ligand_objects)

    for prot_obj in protein_objects:
        # Select interface residues
        interface_sel = f"interface_{prot_obj}"
        cmd.select(
            interface_sel,
            f"(byres ({prot_obj} within {INTERFACE_DIST} of ({all_ligands}))) and {prot_obj}"
        )

        # Show as lines (in addition to cartoon)
        cmd.show("lines", interface_sel)

        # Hide hydrogens on interface residues too
        cmd.hide("everything", f"{interface_sel} and hydro")

        n_res = cmd.count_atoms(f"{interface_sel} and name CA")
        print(f"Interface residues for {prot_obj}: {n_res} residues within {INTERFACE_DIST} A")

        # Disable selection to clean up
        cmd.disable(interface_sel)

def finalize_session():
    """Final adjustments and save session."""
    # Center and zoom
    cmd.center("all")
    cmd.zoom("all", buffer=2)

    # Orient for best view
    cmd.orient("all")

    # Deselect everything
    cmd.deselect()

    # Save session
    cmd.save(OUTPUT_FILE)
    print(f"\nSaved PyMOL session to: {OUTPUT_FILE}")

def main():
    """Main function to create PyMOL session."""
    print("=== PyMOL Session Generator ===\n")

    # Setup visualization settings
    setup_visualization()

    # Load structures
    protein_objects = load_proteins()
    ligand_objects = load_ligands()

    # Align proteins
    if DO_ALIGN and len(protein_objects) > 1:
        print("\n--- Aligning structures ---")
        align_structures(protein_objects)

    # Apply styles
    print("\n--- Applying visualization ---")
    style_proteins(protein_objects)
    style_ligands(ligand_objects)
    show_interface_residues(protein_objects, ligand_objects)

    # Finalize and save
    print("\n--- Finalizing ---")
    finalize_session()

if __name__ == "__main__":
    main()
```

### 4. Run PyMOL

Execute the script with PyMOL:

```bash
# Headless mode (no GUI)
pymol -cq pymol_script.py

# Or with GUI to verify
pymol pymol_script.py
```

### 5. Alternative: Direct Command Generation

If PyMOL is not installed or user prefers commands, output a `.pml` script instead:

```pml
# PyMOL script for protein-ligand visualization
# Generated by /pymol command

# Settings
bg_color black
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1
set stick_radius, 0.15
set line_width, 2

# Load proteins
load protein1.pdb, protein1
load protein2.pdb, protein2

# Load ligands
load ligand.sdf, ligand

# Align structures
align protein2 and name CA, protein1 and name CA

# Style proteins - cartoon
hide everything
show cartoon, protein1
show cartoon, protein2
util.cbc protein1
util.cbc protein2

# Style ligands - sticks, no hydrogens
show sticks, ligand
hide everything, ligand and hydro
util.cbag ligand
set stick_radius, 0.2, ligand

# Interface residues - lines within 5A of ligand
select interface, byres (protein1 within 5.0 of ligand)
show lines, interface
hide everything, interface and hydro
disable interface

# Remove hydrogens from proteins
hide everything, protein1 and hydro
hide everything, protein2 and hydro

# Final view
center all
zoom all, 2
orient all
deselect

# Save session
save session.pse
```

### 6. Report Results

After creating the session, report:
1. Number of protein structures loaded
2. Number of ligand poses loaded
3. Alignment RMSDs (if applicable)
4. Number of interface residues identified
5. Output file path

## Ensemble Handling

When `--ensemble` is specified or multiple PDBs from same directory:

```python
# Detect if files are an ensemble (e.g., backrub output)
if all files have similar names with numbering:
    # Load as states
    cmd.load(pdb_files[0], "ensemble")
    for i, f in enumerate(pdb_files[1:], start=2):
        cmd.load(f, "ensemble", state=i)
```

Navigate states in PyMOL:
- `set all_states, 1` - Show all states superimposed
- `mplay` - Animate through states
- `frame N` - Go to state N

## Color Schemes

| Scheme | Description |
|--------|-------------|
| `chain` | Different color per chain (default) |
| `spectrum` | Rainbow from N to C terminus |
| `ss` | Red=helix, Yellow=sheet, Green=loop |

## Troubleshooting

### PyMOL not found
```bash
# Install via conda
conda install -c conda-forge pymol-open-source

# Or pip (open-source version)
pip install pymol-open-source
```

### Large ensembles slow
- Use `--no-align` to skip alignment
- Load subset: `ensemble_00{01..10}.pdb`

### Ligand not visible
- Check ligand file format (SDF preferred)
- Verify coordinates overlap with protein

### Session file too small (< 10KB)
If the .pse file is only 1-2 KB, the structures failed to load. Common causes:
- **Wrong file paths**: Always use absolute paths in .pml scripts
- **Path typos**: Watch for underscores vs hyphens in directory names
- **Files don't exist**: Verify with `ls -la` before running PyMOL

**Verify session loaded correctly:**
```bash
# Good session size (should be 100KB+ for a protein)
ls -la session.pse
# -rw-r--r--  690835 Jan  8 14:06 session.pse  # Good!
# -rw-r--r--    1158 Jan  8 14:05 session.pse  # Bad - files didn't load
```

### Selection errors in .pml scripts
PyMOL selection syntax quirks:
```pml
# CORRECT - no parentheses needed for simple within
select interface, byres protein within 5.0 of ligand

# ALSO WORKS - with parentheses
select interface, byres (protein within 5.0 of ligand)

# For hiding, use parentheses around compound selections
hide (interface and hydro)
```

### Python scripts (.py) vs PML scripts (.pml)
**Prefer .pml scripts** - they are more reliable:
- `.pml`: Native PyMOL commands, always works
- `.py`: Requires correct Python environment, may have import issues

```bash
# .pml is more reliable
pymol -cq script.pml

# .py may have issues with pymol imports
pymol -cq script.py  # Can fail silently
```

## Best Practices for PyMOL Scripting

### 1. Always use absolute paths
```pml
# GOOD - absolute path
load /Users/name/project/protein.pdb, protein

# BAD - relative path (may fail depending on working directory)
load protein.pdb, protein
```

### 2. Get the working directory path correctly
```bash
# Check actual path (watch for underscores vs hyphens!)
pwd
# /Users/name/project_name  <- underscores
# NOT /Users/name/project-name  <- hyphens
```

### 3. Verify files exist before scripting
```bash
ls -la /full/path/to/files/*.pdb
```

### 4. Standard .pml template
```pml
# Settings
bg_color black
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1
set stick_radius, 0.15
set line_width, 2

# Load with ABSOLUTE paths
load /absolute/path/to/protein.pdb, protein
load /absolute/path/to/ligand.sdf, ligand

# Style protein
hide everything
show cartoon, protein
util.cbc protein
hide (protein and hydro)

# Style ligand
show sticks, ligand
hide (ligand and hydro)
util.cbag ligand
set stick_radius, 0.2, ligand

# Interface residues
select interface, byres protein within 5.0 of ligand
show lines, interface
hide (interface and hydro)
disable interface

# Final view
center all
zoom all, 2
orient all
deselect

# Save with ABSOLUTE path
save /absolute/path/to/session.pse
```

### 5. Running PyMOL

```bash
# Headless mode - create session without GUI
pymol -cq script.pml

# Open session with GUI (macOS)
open -a PyMOL session.pse

# Open session with GUI (Linux/generic)
pymol session.pse
```

### 6. Useful PyMOL commands for inspection
```pml
# Show all loaded objects
get_names

# Count atoms in selection
count_atoms protein
count_atoms ligand

# Check number of states (for multi-pose SDF)
count_states ligand

# Show all states superimposed
set all_states, 1

# Animate through states
mplay
```

## Output

The command generates:
1. `session.pse` - PyMOL session file (or custom name)
2. Console output with alignment RMSDs and interface residue counts

Open the session:
```bash
# macOS - opens PyMOL GUI
open -a PyMOL session.pse

# Linux/generic
pymol session.pse
```
