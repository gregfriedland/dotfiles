---
name: atlas-rank
description: Use when the user asks about "rank benchmark", "ranking benchmark", "LIT-PCBA", "DEKOIS", "virtual screening benchmark", running docking benchmarks, analyzing benchmark results, or evaluating compound scoring workflows.
version: 1.1.0
---

# Atlas Rank Benchmark

Benchmark system for evaluating protein-compound scoring/docking workflows against standard datasets (LIT-PCBA, DEKOIS 2.0).

## Location

`rezo/atlas/benchmarks/rank/` (note: "benchmarks" with an 's')

## Architecture

```
rezo/atlas/benchmarks/rank/
├── benchmark.py          # Core RankBenchmarkConfig class
├── bin/__main__.py       # CLI entry point with 4 commands
├── compound_loader.py    # Loads compounds from GCS
├── constants.py          # Enums, column names, target lists
├── db.py                 # BigQuery database operations
├── flyte.py              # Flyte workflow for parallel execution
├── scorers/
│   ├── base.py           # Abstract CompoundBenchmarkScorer
│   └── dock.py           # RankBenchmarkDockScorer implementation
└── tests/
```

## Supported Benchmarks

### LIT-PCBA (default)
- Data: `gs://rezo-atlas/datasets/LIT_PCBA/`
- 15 targets with crystal structures
- Targets: ADRB2, ALDH1, ESR1_ago, ESR1_ant, FEN1, GBA, IDH1, KAT2A, MAPK1, MTORC1, OPRK1, PKM2, PPARG, TP53, VDR
- Reference: https://pubs.acs.org/doi/10.1021/acs.jcim.0c00155

### DEKOIS 2.0
- Data: `gs://rezo-atlas/datasets/DEKOIS_2.0/`
- 81 targets with challenging physicochemically-matched decoys
- Reference: https://pubs.acs.org/doi/10.1021/ci400115b

## Database Schema

BigQuery dataset: `rezo-atlas.atlas_benchmarks`

| Table | Purpose |
|-------|---------|
| `rank_actives` | Actives/decoys with SMILES and `is_active` flag |
| `rank_scorers` | Scorer configurations with `rank_col` and `rank_col_lower_is_better` |
| `rank_scores` | Compound scores per target |

## Metrics Calculated

- **EF@1%, EF@5%**: Enrichment factors at early percentages
- **AUROC**: Area Under ROC Curve
- **BEDROC_20**: BEDROC with alpha=20 (early enrichment focus)
- **AP**: Average Precision
- **Log-AUC**: Area under semi-log ROC curve

## CLI Commands

### 1. Run a Benchmark

#### Basic usage with Rosetta docking:

```bash
python -m rezo.atlas.benchmarks.rank.bin run \
    --scorer.rosetta-dock.mode vsh \
    --scorer.rosetta-dock.num-poses 3 \
    --scorer.rosetta-dock.user greg \
    --rank-col gnina_combined \
    --targets ADRB2 ALDH1
```

#### Using GNINA docking:

```bash
python -m rezo.atlas.benchmarks.rank.bin run \
    --scorer.gnina-dock.exhaustiveness 1000 \
    --scorer.gnina-dock.num-modes 20 \
    --scorer.gnina-dock.num-mc-steps 100 \
    --scorer.gnina-dock.box-size 20 \
    --scorer.gnina-dock.cnn-scoring rescore \
    --scorer.gnina-dock.cnn-model default \
    --scorer.gnina-dock.cnn-scoring-mult 20 \
    --scorer.gnina-dock.seed 42 \
    --scorer.gnina-dock.use-gpu false \
    --scorer.score-cols vina_affinity gnina_pose_score \
    --targets ADRB2 ALDH1
```

#### Key options:

**Benchmark options:**
- `--benchmark`: `litpcba` (default) or `dekois2`
- `--targets`: Specific targets (default: all targets for benchmark)
- `--max-decoys-per-active`: Max decoys per active (default: 50)
- `--random-seed`: For reproducibility (default: 42)

**Ranking options (part of scorer config):**
- `--scorer.rank-col`: Score column for ranking (default: `aggregate_rank`)
  - `aggregate_rank`: Computed from all score_cols during pose aggregation (recommended)
  - For GNINA docking: `vina_affinity`, `gnina_pose_score`, `gnina_affinity`, `gnina_combined`
  - For Rosetta docking: `docking_score`

**GNINA docking options** (`--scorer.gnina-dock.*`):
- `--scorer.gnina-dock.exhaustiveness`: Search exhaustiveness (default: 1000)
- `--scorer.gnina-dock.num-modes`: Number of binding modes (default: 5)
- `--scorer.gnina-dock.num-mc-steps`: Monte Carlo steps per ligand
- `--scorer.gnina-dock.box-size`: Docking box size in Angstroms (default: 20)
- `--scorer.gnina-dock.cnn-scoring`: CNN scoring mode (`none`, `rescore`, `refinement`, `metrorescore`, `metrorefine`, `all`)
- `--scorer.gnina-dock.cnn-model`: CNN model (`default`, `fast`, `dense`, etc.)
- `--scorer.gnina-dock.cnn-scoring-mult`: CNN rescoring multiplier (default: 20)
- `--scorer.gnina-dock.seed`: Random seed (default: 42)
- `--scorer.gnina-dock.use-gpu`: Use GPU (default: true)
- `--scorer.gnina-dock.fast-line-search`: Fast GPU BFGS line search (default: false)
- `--scorer.gnina-dock.mc-bfgs-iters-hunt`: BFGS iterations for MC hunt phase
- `--scorer.gnina-dock.mc-bfgs-iters-refine`: BFGS iterations for MC refine phase
- `--scorer.gnina-dock.post-mc-bfgs-iters-refine`: Post-MC BFGS iterations

**Rosetta docking options** (`--scorer.rosetta-dock.*`):
- `--scorer.rosetta-dock.mode`: Docking mode (`vsh`, `random`, etc.)
- `--scorer.rosetta-dock.num-poses`: Number of poses (default: 3)
- `--scorer.rosetta-dock.user`: User for docking job

**Other options:**
- `--medchem-alert-filter`: Filter compounds with alerts (PAINS, BMS, Glaxo)
- `--force-rerun`: Bypass deduplication and create new scorer

### 2. Analyze Results

```bash
python -m rezo.atlas.benchmarks.rank.bin analyze \
    --scorer-ids scorer_abc123 scorer_def456 \
    --benchmark litpcba \
    --pairwise \
    --min-seeds 3
```

Options:
- `--scorer-ids`: One or more scorer IDs to analyze
- `--benchmark`: `litpcba` or `dekois2`
- `--pairwise`: Run Wilcoxon signed-rank tests between config groups
- `--min-seeds`: Minimum seeds required per config for pairwise comparison
- `--verbose`: Show per-target results

Output:
- Summary table with median logAUC and BEDROC per scorer
- Per-target results (with `--verbose`)
- Pairwise comparisons showing statistical significance

### 3. List All Scorers

```bash
python -m rezo.atlas.benchmarks.rank.bin list \
    --benchmark litpcba \
    --target ADRB2 \
    --limit 100
```

Options:
- `--benchmark`: `litpcba` or `dekois2`
- `--target`: Filter to scorers with data for specific target
- `--limit`: Maximum scorers to display (default: 100)

Shows: scorer IDs, configs, progress, error counts, completion status.

### 4. Populate Actives Table (One-time Setup)

```bash
python -m rezo.atlas.benchmarks.rank.bin load --benchmark litpcba
```

Loads all compound metadata from the dataset into BigQuery.

## Execution Flow

1. **CLI** parses args using `PydanticArgparse`
2. **`run_command()`** creates:
   - `RankBenchmarkDB` for database operations
   - `RankBenchmarkConfig` with benchmark settings
   - `RankBenchmarkDockScorer` from scorer config + rank_col
3. **Flyte execution**: `run_benchmark()` workflow runs on `ATLAS_WEST1`
   - Creates/reuses scorer_id based on config hash
   - Runs `_run_single_target_benchmark()` in parallel for each target
   - Each target: loads compounds -> scores -> inserts to BigQuery
4. **Analysis**: Fetches scores from BigQuery, calculates metrics

## Key Classes

| Class | Purpose |
|-------|---------|
| `RankBenchmarkConfig` | Main benchmark config with target validation |
| `RankBenchmarkDB` | BigQuery CRUD operations |
| `RankBenchmarkCompoundLoader` | Loads compounds from GCS with sampling |
| `CompoundBenchmarkScorer` | Abstract base for scorers |
| `RankBenchmarkDockScorer` | Docking scorer using `score_compounds_dynamic` |

## Available Score Columns (rank_col)

**Default (recommended):**
- `aggregate_rank`: Combined rank from all score_cols (lower is better)

For **GNINA docking** (`--scorer.gnina-dock.*`):
- `vina_affinity`: Vina binding affinity (kcal/mol, lower is better)
- `gnina_pose_score`: CNN pose score (higher is better)
- `gnina_affinity`: CNN affinity prediction (lower is better)
- `gnina_combined`: Combined CNN affinity + pose score (higher is better)

For **Rosetta docking** (`--scorer.rosetta-dock.*`):
- `docking_score`: Rosetta docking score (lower is better)

For **Rosetta + GNINA rescore** (`--scorer.gnina-rescore.*`):
- All GNINA scores above, plus `docking_gnina_combined`

## Baseline Settings

The analyze command compares scorers against this baseline:

```python
BASELINE_SETTINGS = {
    "mode": "vsh",
    "num_poses": 3,
    "partial_charge": "antechamber",
    "relax": False,
    "hybridize": False,
    "score_col": "gnina_combined",
    "pose_agg": "mean",
}
```

Scorers matching these settings are labeled "(baseline)" in analysis output.

## Example Workflows

### Run benchmark with GNINA docking (CPU, CNN rescore)

```bash
python -m rezo.atlas.benchmarks.rank.bin run \
    --scorer.gnina-dock.exhaustiveness 1000 \
    --scorer.gnina-dock.num-modes 20 \
    --scorer.gnina-dock.cnn-scoring rescore \
    --scorer.gnina-dock.cnn-model default \
    --scorer.gnina-dock.use-gpu false \
    --targets ADRB2 ALDH1
```

### Run benchmark with GNINA docking (pure Vina scoring)

```bash
python -m rezo.atlas.benchmarks.rank.bin run \
    --scorer.gnina-dock.exhaustiveness 1000 \
    --scorer.gnina-dock.num-modes 10 \
    --scorer.gnina-dock.cnn-scoring none \
    --targets ADRB2
```

### Run benchmark on all LIT-PCBA targets with Rosetta

```bash
python -m rezo.atlas.benchmarks.rank.bin run \
    --scorer.rosetta-dock.mode vsh \
    --scorer.rosetta-dock.num-poses 3 \
    --scorer.rosetta-dock.user greg
```

### Compare multiple scorers

```bash
python -m rezo.atlas.benchmarks.rank.bin analyze \
    --scorer-ids scorer_abc123 scorer_def456 scorer_ghi789 \
    --pairwise \
    --verbose
```

## Troubleshooting

### Finding scorer IDs
Use `list` command to see all scorers and their IDs.

### Checking progress
The `list` command shows how many targets have been completed for each scorer.

### Resuming failed runs
The system automatically skips targets that already have scores. Just re-run with the same config to resume.

### Force rerun
Use `--force-rerun` to create a new scorer even if identical config exists.

### Invalid rank_col error
Make sure the `rank_col` is available from your docking configuration:
- GNINA docking provides: `vina_affinity`, `gnina_pose_score`, `gnina_affinity`, `gnina_combined`
- Rosetta docking provides: `docking_score`
- Adding `--scorer.gnina-rescore.*` to Rosetta adds GNINA scores
