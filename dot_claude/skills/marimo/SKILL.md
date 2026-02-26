---
name: marimo
description: Use when the user asks to create a marimo notebook script (.mo.py). Examples - "create a marimo notebook", "new marimo notebook", "make a .mo.py file".
---

# Marimo Notebook Creator

When the user asks you to create a marimo notebook, generate a `.mo.py` file following the exact patterns below.

## Arguments

The user invokes `/marimo <filename>`. If `$ARGUMENTS` is provided, use it as the filename. The file should end in `.mo.py`. Create it in the current working directory unless the user specifies otherwise.

## Core Structure

Every marimo notebook MUST start with this exact boilerplate:

```python
import marimo

__generated_with = "0.19.11"
app = marimo.App(width="medium")
```

And end with:

```python
if __name__ == "__main__":
    app.run()
```

## Cell Rules

- Each cell is a function decorated with `@app.cell`.
- **ALL imports MUST go in the first cell** (named `def imports():`). This avoids marimo variable redefinition errors.
- Cells return a tuple of the variables they export (e.g., `return (compounds, df)`). If nothing is exported, use `return` with no value.
- Cells that use variables from other cells receive them as function arguments.
- Prefix temporary/local variables with `_` to avoid marimo tracking them as cell outputs.
- Use `mo.md(...)` for markdown output, `mo.vstack([...])` to stack elements vertically.
- Use `mo.ui.table(df)` to display interactive tables.

## First Cell: Imports + Disk Cache + Constants

The first cell MUST always be named `def imports():` and contain:

1. All standard library imports
2. All third-party imports
3. The rezo monorepo sys.path setup
4. All `from rezo.*` imports needed by the notebook
5. The `_make_hashable` and `disk_cache` utility functions
6. Constants (NOTEBOOK_DIR, file paths, etc.)

Here is the template for the first cell:

```python
@app.cell
def imports():
    # All imports should go in this cell to avoid marimo variable redefinition errors.
    import functools
    import hashlib
    import math
    import os
    import pickle
    import re
    import subprocess
    import sys
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd
    import plotly.express as px
    from rdkit import Chem

    # --- Rezo monorepo setup ---
    repo_root = Path("/Users/gregfriedland/src/rezo")
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from rezo.common.biochem.compounds.cluster import (
        MolClusterConfig,
        MolClustererType,
    )
    from rezo.common.biochem.fingerprint import FPSimilarityConfig
    from rezo.common.biochem.compounds.compounds1d import (
        Compounds1D,
    )
    from rezo.common.biochem.compounds.properties.descriptor import (
        RdkitDescriptorProperties,
    )
    from rezo.common.biochem.compounds.properties.medchem_alert import (
        MEDCHEM_ALERT_PREFIX,
        MedchemAlertProperties,
    )
    from rezo.common.biochem.draw import draw_image_grid_with_alignment

    # --- Disk cache utilities ---
    def _make_hashable(obj):
        """Convert an object to bytes for cache key hashing."""
        if isinstance(obj, pd.DataFrame):
            return pd.util.hash_pandas_object(obj).values.tobytes()
        if hasattr(obj, "df") and isinstance(getattr(obj, "df"), pd.DataFrame):
            return pd.util.hash_pandas_object(obj.df).values.tobytes()
        try:
            return pickle.dumps(obj)
        except Exception:
            return repr(obj).encode()

    def disk_cache(cache_dir, name=None):
        """Decorator that caches function results to disk based on hashed inputs.

        Usage::

            @disk_cache(NOTEBOOK_DIR, name="descriptors")
            def calc(compounds):
                ...
                return result

            result = calc(compounds)
            # calc.was_cached is True/False, calc.cache_path is the file path
        """
        def decorator(fn):
            _name = name or fn.__name__

            @functools.wraps(fn)
            def wrapper(*args, **kwargs):
                h = hashlib.sha256()
                h.update(_name.encode())
                for arg in args:
                    h.update(_make_hashable(arg))
                for k in sorted(kwargs):
                    h.update(k.encode())
                    h.update(_make_hashable(kwargs[k]))
                cache_key = h.hexdigest()[:12]
                cache_path = cache_dir / f".cache_{_name}_{cache_key}.pkl"

                if cache_path.exists():
                    with open(cache_path, "rb") as f:
                        wrapper.was_cached = True
                        wrapper.cache_path = cache_path
                        return pickle.load(f)

                result = fn(*args, **kwargs)
                with open(cache_path, "wb") as f:
                    pickle.dump(result, f)
                wrapper.was_cached = False
                wrapper.cache_path = cache_path
                return result

            wrapper.was_cached = False
            wrapper.cache_path = None
            return wrapper
        return decorator

    # --- Constants ---
    NOTEBOOK_DIR = Path(__file__).parent
    VENV_PYTHON = str(repo_root / ".venv" / "bin" / "python")
    ANSI_ESCAPE = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")

    # Return all symbols that other cells need
    return (
        Chem,
        Compounds1D,
        FPSimilarityConfig,
        MEDCHEM_ALERT_PREFIX,
        MedchemAlertProperties,
        MolClusterConfig,
        MolClustererType,
        NOTEBOOK_DIR,
        RdkitDescriptorProperties,
        disk_cache,
        draw_image_grid_with_alignment,
        mo,
        np,
        pd,
        px,
        repo_root,
    )
```

**IMPORTANT**: Adjust the `repo_root` path to point to whichever rezo worktree or checkout the user is working with. Ask if unclear. Also adjust the return tuple to only include symbols actually needed by the notebook's cells.

## Second Cell: Title Markdown

```python
@app.cell
def _(mo):
    mo.md("""
    # Notebook Title Here

    Brief description of what this notebook does.
    """)
    return
```

## Pattern: Molecular Descriptors

Calculate RDKit molecular descriptors on a `Compounds1D` object:

```python
@app.cell
def calc_descriptors(
    NOTEBOOK_DIR,
    RdkitDescriptorProperties,
    compounds_std,
    disk_cache,
    mo,
):
    """Calculate RDKit molecular descriptors, caching to disk."""
    DESCRIPTOR_PROPS = ["MolWt", "NumRotatableBonds", "MolLogP", "QED"]

    @disk_cache(NOTEBOOK_DIR, name="descriptors")
    def _calc(compounds):
        desc_calc = RdkitDescriptorProperties(properties=DESCRIPTOR_PROPS)
        return desc_calc.calculate(compounds)

    compounds_desc = _calc(compounds_std)
    _status = "cached" if _calc.was_cached else "calculated"
    mo.md(f"Descriptors **{_status}**: {', '.join(DESCRIPTOR_PROPS)} (`{_calc.cache_path.name}`)")
    return DESCRIPTOR_PROPS, compounds_desc
```

## Pattern: Medchem Alerts

Run medchem alert filters on a `Compounds1D` object:

```python
@app.cell
def calc_alerts(
    MEDCHEM_ALERT_PREFIX,
    MedchemAlertProperties,
    NOTEBOOK_DIR,
    compounds_desc,
    disk_cache,
    mo,
):
    """Run medchem alert filters, caching to disk."""
    ALERT_SETS = ["PAINS", "Dundee", "BMS", "Reactive-Unstable-Toxic"]
    ALERT_COLS = [f"{MEDCHEM_ALERT_PREFIX}{a}" for a in ALERT_SETS]

    @disk_cache(NOTEBOOK_DIR, name="alerts")
    def _calc(compounds_desc):
        alert_calc = MedchemAlertProperties(properties=ALERT_COLS)
        return alert_calc.calculate(compounds_desc)

    compounds_scored = _calc(compounds_desc)
    _status = "cached" if _calc.was_cached else "calculated"
    mo.md(f"Alerts **{_status}**: {', '.join(ALERT_SETS)} (`{_calc.cache_path.name}`)")
    return ALERT_COLS, ALERT_SETS, compounds_scored
```

## Pattern: Custom Enumeration

Run custom library enumeration using `enumerate_custom_library`. This is typically done via a separate script in a `scripts/` subdirectory, invoked from a notebook cell. The script pattern:

```python
# scripts/enumerate_XYZ.py
from rezo.atlas.library.custom.enumerate import enumerate_custom_library
from rezo.common.biochem.compounds.compounds1d import Compounds1D

products = enumerate_custom_library(
    target_smiles=INTERMEDIATE_SMILES,      # SMILES of the scaffold with reactive handle
    target_id="INTERMEDIATE",                # identifier for the scaffold
    building_blocks_paths=["enamine_stock", "wuxi_parallel"],  # BB sources
    config_name="buchwald",                  # reaction config name (e.g., "buchwald", "suzuki")
    reaction_ids=["cn_primary_arbr", ...],   # specific reaction IDs to use
    reverse_components=False,                # whether to swap target/BB roles
    max_products=None,                       # optional cap on products
)
products.to_file("products.parquet")
```

The notebook cell to invoke enumeration (typically commented out after first run):

```python
@app.cell
def run_enumeration():
    """Run enumeration script (currently disabled).

    To run manually:
        cd <repo_root>
        PYTHONPATH=. VENV_PATH=.venv python scripts/enumerate_XYZ.py \\
            --output_file <PRODUCTS_FILE>
    """
    # import os
    # result = subprocess.run(
    #     [VENV_PYTHON, str(NOTEBOOK_DIR / "scripts" / "enumerate_XYZ.py"),
    #      "--output_file", str(PRODUCTS_FILE)],
    #     capture_output=True, text=True,
    #     cwd=str(repo_root),
    #     env={**os.environ, "PYTHONPATH": str(repo_root)},
    #     timeout=600,
    # )
    return
```

## Pattern: GNINA Docking

Submit GNINA docking via the ATLAS API with disk caching:

```python
@app.cell
def run_docking(
    NOTEBOOK_DIR,
    RECEPTOR_PDB,
    compounds_filtered,
    disk_cache,
    mo,
):
    """Submit GNINA docking via Flyte, with disk caching."""
    from rezo.atlas.api import ATLASApi
    from rezo.common.biochem.target3d import BindingSiteCenter, Target3D
    from rezo.atlas.library.explorer.base.scorer_config import CompoundScorerConfig
    from rezo.common.biochem.compounds.compounds1d import TargetCompounds1D
    from rezo.atlas_models.gnina.config import (
        GninaCNNModel,
        GninaCNNScoring,
        GninaDockConfig,
        GninaRefinementMode,
        GninaRescoringConfig,
    )

    gnina_config = GninaDockConfig(
        exhaustiveness=10_000,
        num_mc_steps=100,
        num_modes=10,
        use_gpu=False,
        cnn_scoring=GninaCNNScoring.RESCORE,
        cnn_model=GninaCNNModel.DEFAULT,
    )
    gnina_rescore_config = GninaRescoringConfig(
        refinement_mode=GninaRefinementMode.MINIMIZE,
        cnn_scoring=GninaCNNScoring.REFINEMENT,
    )
    scorer_config = CompoundScorerConfig(
        gnina_dock=gnina_config,
        gnina_rescore=gnina_rescore_config,
    )

    DOCKED_FILE = NOTEBOOK_DIR / "docked.parquet"

    @disk_cache(NOTEBOOK_DIR, name="gnina_docking")
    def _dock(compounds):
        target = Target3D.from_file(RECEPTOR_PDB, name="TARGET")
        target.binding_site_center = BindingSiteCenter.from_file(
            RECEPTOR_PDB.with_suffix(".bsxyz")
        )
        target_compounds = TargetCompounds1D.from_tc_pair(
            target, compounds
        )
        return ATLASApi.score(
            config=scorer_config,
            compounds=target_compounds,
            wait=True,
        )

    docked = _dock(compounds_filtered)
    docked.to_file(DOCKED_FILE)
    _status = "cached" if _dock.was_cached else "submitted"
    mo.md(
        f"Docking **{_status}**: **{len(docked)}** poses "
        f"(exh={gnina_config.exhaustiveness}, nmc={gnina_config.num_mc_steps}) "
        f"saved to `{DOCKED_FILE.name}`"
    )
    return (docked,)
```

**Notes for docking cells:**
- The receptor PDB and `.bsxyz` file must exist. The `.bsxyz` is auto-found as `RECEPTOR_PDB.with_suffix(".bsxyz")`.
- `ATLASApi.score(wait=True)` blocks until the Flyte workflow finishes.
- Adjust `exhaustiveness`, `num_mc_steps`, and `num_modes` as needed.
- The `FlyteProject` enum can be passed to `ATLASApi.score(project=...)` if needed.

## Pattern: Butina Clustering

Cluster compounds by Tanimoto fingerprint similarity:

```python
@app.cell
def cluster_compounds(
    Compounds1D,
    FPSimilarityConfig,
    MolClusterConfig,
    MolClustererType,
    docked,
    mo,
):
    """Butina-cluster compounds at a given Tanimoto similarity threshold."""
    _cluster_config = MolClusterConfig(
        type_=MolClustererType.BUTINA,
        fp_config=FPSimilarityConfig(similarity_threshold=0.8),
    )
    # Sort by score before clustering so cluster exemplars are best-scoring
    _df = docked.df.sort_values("gnina_affinity", ascending=False)
    _c = Compounds1D.from_df(_df)
    clustered = _cluster_config.cluster(_c)
    _n_clusters = clustered.df["cluster_id"].nunique()
    mo.md(
        f"Clustered **{len(clustered)}** compounds into "
        f"**{_n_clusters}** clusters at 0.8 Tanimoto similarity"
    )
    return (clustered,)
```

## Pattern: Load Compounds from Parquet

```python
@app.cell
def load_compounds(Compounds1D, mo):
    """Load compounds from parquet file."""
    COMPOUNDS_FILE = NOTEBOOK_DIR / "compounds.parquet"
    compounds = Compounds1D.from_file(COMPOUNDS_FILE)
    mo.md(f"Loaded **{len(compounds)}** compounds from `{COMPOUNDS_FILE.name}`")
    return (compounds,)
```

## Pattern: Standardize & Deduplicate

```python
@app.cell
def _(Compounds1D, NOTEBOOK_DIR, compounds, disk_cache, mo):
    """Standardize SMILES and deduplicate."""

    @disk_cache(NOTEBOOK_DIR, name="standardize")
    def _standardize(compounds):
        c = Compounds1D.from_df(compounds.df.copy())
        c.standardize()
        c.deduplicate()
        return c

    _n_before = len(compounds)
    compounds_std = _standardize(compounds)
    _n_after = len(compounds_std)
    _status = "cached" if _standardize.was_cached else "computed"
    mo.md(
        f"Standardization **{_status}**: **{_n_before}** → **{_n_after}** compounds "
        f"(`{_standardize.cache_path.name}`)"
    )
    return (compounds_std,)
```

## Pattern: Score Distribution Plots with Reference Line

```python
@app.cell
def _(go):
    """Define reusable plot_score_distributions helper."""
    import plotly.graph_objects as go

    def plot_score_distributions(df, ref_id, score_cols, title_prefix="Distribution", label="compounds"):
        _ref_mask = df["compound_id"] == ref_id
        _ref_row = df[_ref_mask].iloc[0] if _ref_mask.any() else None
        _lib_df = df[~_ref_mask]
        _figs, _summaries = [], []
        for _col in score_cols:
            _fig = go.Figure()
            _fig.add_trace(go.Histogram(x=_lib_df[_col].dropna(), nbinsx=80, name=label.capitalize(), marker_color="steelblue", opacity=0.8))
            if _ref_row is not None:
                _ref_val = _ref_row[_col]
                _fig.add_vline(x=_ref_val, line_dash="dash", line_color="red", line_width=2,
                    annotation_text=f"{ref_id}: {_ref_val:.3f}", annotation_position="top right")
                _valid = _lib_df[_col].dropna()
                _n_better = int((_valid > _ref_val).sum())
                _pct = 100 * _n_better / len(_valid) if len(_valid) > 0 else 0
                _summaries.append(f"**{_col}**: {_n_better}/{len(_valid)} ({_pct:.1f}%) {label} score higher than {ref_id} ({_ref_val:.3f})")
            _fig.update_layout(title=f"{title_prefix} — {_col}", xaxis_title=_col, yaxis_title="Count", bargap=0.05, height=400)
            _figs.append(_fig)
        return _figs, _summaries
    return (plot_score_distributions,)
```

## Key Reminders

- **Never put imports in cells other than the first `imports()` cell.** The only exception is deferred imports inside functions (e.g., `from rezo.atlas.api import ATLASApi` inside a docking cell function body is OK since it's a heavy import not always needed).
- **Always use `disk_cache`** for any computation that takes more than a few seconds.
- **Always report cache status** in the cell's `mo.md(...)` output so the user knows if results were cached or freshly computed.
- **Use `_` prefix** for temporary variables to avoid marimo accidentally exporting them.
- **Section headers** should be their own `mo.md(...)` cells for clean notebook structure.
