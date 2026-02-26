---
name: atlas-pose
description: Run and analyze pose prediction benchmarks (GNINA, Asterix docking) on posebusters or runs_n_poses datasets via Flyte.
version: 1.0.0
---

# Atlas Pose Benchmark

Run pose prediction benchmarks on posebusters or runs_n_poses datasets using GNINA or Asterix docking.

## Location

The benchmark CLI is in the rezo worktree:
```
~/src/worktrees/greg-260101-pose-benchmark3/
```

## Commands

### Run Benchmark

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin run [OPTIONS]
```

### Analyze Results

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin analyze --benchmark posebusters
```

### Deactivate Scorers

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin deactivate --scorer-ids ID1 ID2
```

### Analyze Specific Scorers

Use short IDs without `pose_scorer_` prefix:

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin analyze \
  --benchmark posebusters \
  --include-ids ID1 ID2 ID3
```

### Get Full Scorer Config (Python)

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -c "
from pprint import pprint
from rezo.atlas.benchmarks.pose.db import PoseBenchmarkDB
db = PoseBenchmarkDB()
config = db.get_scorer_config('pose_scorer_SCORER_ID')
pprint(config)
"
```

## Run Command Flags

### General Options

| Flag | Description | Default |
|------|-------------|---------|
| `--dataset` | Benchmark dataset (`posebusters` or `runs_n_poses`) | `posebusters` |
| `--seed` | Random seed for reproducibility | `42` |
| `--wait` / `--no-wait` | Wait for Flyte execution to complete | `--wait` |
| `--verbose` / `--no-verbose` | Enable verbose logging | `--no-verbose` |
| `--force-rerun` / `--no-force-rerun` | Force rerun even if results exist | `--no-force-rerun` |
| `--max-systems` | Maximum number of systems to evaluate | None (all) |
| `--similarity-bins` | Filter to specific similarity bins (runs_n_poses only) | None |
| `--flyte-project` | Flyte project to submit to | `ATLAS_CENTRAL1` |
| `--scorer.score-cols` | Score columns to use for ranking | Required |
| `--scorer.keep-pose-columns` / `--no-scorer.keep-pose-columns` | Keep pose columns in output | `False` |

### Multi-Value Flags

For flags that accept multiple values (like `--scorer.score-cols`), pass all values after a single flag separated by spaces:

```bash
# Correct: single flag with multiple values
--scorer.score-cols vina_affinity gnina_pose_score

# Wrong: do NOT repeat the flag
--scorer.score-cols vina_affinity --scorer.score-cols gnina_pose_score
```

Available score columns:
- `vina_affinity` - Vina/AutoDock scoring function
- `gnina_pose_score` - GNINA CNN pose score
- `gnina_affinity` - GNINA CNN affinity prediction
- `gnina_combined` - Combined GNINA score
- `docking_score` - Rosetta docking score (when using rosetta_dock)

---

## Asterix Dock Config (`--scorer.asterix-dock.*`)

Asterix is the GPU-accelerated fork of GNINA with Adam optimizer-based docking.

| Flag | Description | Default |
|------|-------------|---------|
| `--scorer.asterix-dock.exhaustiveness` | Search exhaustiveness | `1000` |
| `--scorer.asterix-dock.num-modes` | Number of binding modes to generate | `1` |
| `--scorer.asterix-dock.adam-iters` | Adam optimizer iterations per MC step | `100` |
| `--scorer.asterix-dock.refine-iters` | Refinement iterations after MC | `50` |
| `--scorer.asterix-dock.adam-lr` | Adam learning rate | `0.05` |
| `--scorer.asterix-dock.adam-max-step` | Max Adam step size | `1.0` |
| `--scorer.asterix-dock.refine-mult` | Refinement multiplier | `100` |
| `--scorer.asterix-dock.gpu-batch-size` | GPU batch size | `50000` |
| `--scorer.asterix-dock.box-size` | Docking box size in Angstroms | `20.0` |
| `--scorer.asterix-dock.seed` | Random seed | `42` |
| `--scorer.asterix-dock.cnn-scoring` | CNN scoring mode | `none` |
| `--scorer.asterix-dock.cnn-model` | CNN model architecture | None |
| `--scorer.asterix-dock.cnn-scoring-mult` | Poses to rescore with CNN | `20` |
| `--scorer.asterix-dock.max-valid-score` | Max valid score threshold | `1e10` |
| `--scorer.asterix-dock.fast-line-search-trials` | Fast line search trials | `10` |
| `--scorer.asterix-dock.grid-granularity` | Grid granularity | `0.375` |
| `--scorer.asterix-dock.intramol-period` | Intramolecular evaluation period | `5` |
| `--scorer.asterix-dock.keep-pose-columns` / `--no-scorer.asterix-dock.keep-pose-columns` | Keep pose columns | `False` |

---

## GNINA Dock Config (`--scorer.gnina-dock.*`)

GNINA original (CPU-based) docking with Monte Carlo search.

| Flag | Description | Default |
|------|-------------|---------|
| `--scorer.gnina-dock.exhaustiveness` | Search exhaustiveness | `8` |
| `--scorer.gnina-dock.num-modes` | Number of binding modes to generate | `9` |
| `--scorer.gnina-dock.num-mc-steps` | Monte Carlo steps per ligand | None (auto) |
| `--scorer.gnina-dock.box-size` | Docking box size in Angstroms | `20.0` |
| `--scorer.gnina-dock.seed` | Random seed | `42` |
| `--scorer.gnina-dock.cnn-scoring` | CNN scoring mode | `none` |
| `--scorer.gnina-dock.cnn-model` | CNN model architecture | None |
| `--scorer.gnina-dock.keep-pose-columns` / `--no-scorer.gnina-dock.keep-pose-columns` | Keep pose columns | `False` |

### CNN Scoring Modes (`cnn-scoring`)

| Value | Description |
|-------|-------------|
| `none` | No CNN scoring (pure Vina/SMINA) |
| `rescore` | Dock with Vina, rescore final poses with CNN |
| `refinement` | CNN gradient-based minimization after docking |
| `metrorescore` | CNN for Metropolis accept/reject + final scoring |
| `metrorefine` | CNN for Metropolis accept/reject + local refinement |
| `all` | CNN fully replaces Vina during entire search |

### CNN Models (`cnn-model`)

| Value | Description |
|-------|-------------|
| `default` | GNINA's default ensemble (3 models) |
| `fast` | Single fast model for high-throughput |
| `default1.0` | GNINA 1.0 ensemble (5 models, slower) |
| `dense` | Dense model architecture ensemble |
| `dense_1_3` | Single dense model |
| `dense_1_3_PT_KD_3` | Knowledge-distilled dense model |
| `crossdock_default2018` | CrossDocked dataset ensemble |
| `crossdock_default2018_KD_4` | Knowledge-distilled CrossDocked model |
| `general_default2018` | General default 2018 ensemble |
| `redock_default2018` | Redocking data ensemble |

**Note:** `cnn-model` must be set when `cnn-scoring` is not `none`.

---

## Analyze Command Flags

| Flag | Description | Default |
|------|-------------|---------|
| `--benchmark` | Benchmark type (`posebusters` or `runs_n_poses`) | Required |
| `--include-ids` | Specific scorer IDs to include (partial match) | None |
| `--min-coverage` | Minimum coverage percentage (0-100) | `0.0` |
| `--min-systems` | Minimum number of systems | `0` |
| `--exclude` | Scorer IDs to exclude (partial match) | None |
| `--limit` | Maximum scorers to list | `50` |
| `--verbose` / `--no-verbose` | Enable verbose logging | `--no-verbose` |

---

## Examples

### Asterix Docking (GPU, Default Config)

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin run \
  --dataset posebusters \
  --no-wait \
  --scorer.asterix-dock.exhaustiveness 1000 \
  --scorer.asterix-dock.adam-iters 100 \
  --scorer.asterix-dock.refine-iters 50 \
  --scorer.asterix-dock.num-modes 1 \
  --scorer.score-cols vina_affinity
```

### GNINA GPU Docking with CNN Rescoring

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin run \
  --dataset posebusters \
  --no-wait \
  --scorer.gnina-dock.exhaustiveness 1000 \
  --scorer.gnina-dock.num-mc-steps 100 \
  --scorer.gnina-dock.num-modes 20 \
  --scorer.gnina-dock.cnn-scoring rescore \
  --scorer.gnina-dock.cnn-model default \
  --scorer.score-cols vina_affinity gnina_pose_score
```

### GNINA CPU Docking (No CNN)

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin run \
  --dataset posebusters \
  --no-wait \
  --scorer.gnina-dock.exhaustiveness 100 \
  --scorer.gnina-dock.num-mc-steps 100 \
  --scorer.gnina-dock.num-modes 20 \
  --scorer.gnina-dock.cnn-scoring none \
  --scorer.score-cols vina_affinity
```

### Analyze Results

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin analyze \
  --benchmark posebusters \
  --min-coverage 90 \
  --limit 100
```

### Compare Specific Scorers

```bash
cd ~/src/worktrees/greg-260101-pose-benchmark3
PYTHONPATH=. .venv/bin/python -m rezo.atlas.benchmarks.pose.bin analyze \
  --benchmark posebusters \
  --include-ids d94e56765491 8c224cc6d455
```

---

## Output

The `run` command submits a Flyte workflow and returns:
- Execution URL for monitoring
- Scorer ID for analyzing results

The `analyze` command shows a table with:
- `name`: Config summary (non-default values)
- `num_systems`: Total systems evaluated
- `coverage_pct`: Percentage with RMSD computed
- `success_pct`: Percentage with RMSD < 2.0 Angstrom
- `median_rmsd`: Median RMSD across all systems
- `success_pct_no_cof`: Success rate excluding cofactor systems
- `median_rmsd_no_cof`: Median RMSD excluding cofactor systems
- `scorer_id`: Unique identifier
- `created_at`: Timestamp (Pacific time)
