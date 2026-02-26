# Wrapper for running CLI commands

**ALWAYS** use the following wrapper for running CLI commands:
```bash
~/.claude/run_cli.py run "<command> <arguments>"
```
Example:
```bash
~/.claude/run_cli.py run "python3.10 -m infra.environment.bin update --env r2d2"
```


* When using this wrapper, **NEVER** redirect stdout or stderr to a file, since the command does that automatically. This log file is referenced in the output of the command.
* The last 20 lines of the log file are also included in the output of the command.


# How to run python commands

**ALWAYS** run python/pytest commands using the following technique:

* Use the wrapper above to run the command.
* Set the PYTHONPATH to the root of the monorepo and set VENV_PATH to the virtual environment.
* Example:
```bash
~/.claude/run_cli.py run "PYTHONPATH=. VENV_PATH=.venv pytest tests/"
```


# Memory System

You have access to persistent memory via the MCP memory service.

## Protocol
* **On every user request**: call `retrieve_memory` or `search_by_tag`
  with the user's topic to check for relevant prior context.
* **Before responding**: call `store_memory` with any key decisions,
  results, or user preferences learned during this turn.

## What to Store
- Key decisions and outcomes
- User preferences
- Important context for future sessions

## What NOT to Store
- Raw data dumps
- Secrets / API keys
- Temporary debugging info


# Writing Bug Repros

When creating a minimal reproducible example for a bug:

1. **Strip all internal dependencies.** Use only the library under test plus standard/common packages (pandas, pydantic, etc.). No internal modules, no custom wrappers, no project-specific utilities.
2. **One file, copy-pasteable.** The repro should be a single self-contained script that anyone with the library installed can run.
3. **Include a flag to toggle working vs broken behavior.** This makes it trivial to confirm the bug and rules out environmental issues. For example, `--bare` for the working path vs default for the broken path.
4. **Add a timeout with stack dump.** For deadlocks/hangs, use `signal.SIGALRM` to dump all thread stacks after a timeout so the hang location is obvious.
5. **Test both local and remote paths.** Bugs may only manifest in one mode. For Flyte: `flyte run --local` (CLI local), `flyte run` (CLI remote), and `flyte.run()` (Python API remote) are three distinct code paths.
6. **Minimize data.** Use the smallest possible input (e.g., 3-row DataFrame) that still triggers the issue.
7. **Name the file descriptively.** e.g., `test_pydantic_flyte_df_pyapi.py` — include what's being tested and how it's invoked.


# Build Cache Policy

**NEVER DELETE THE BUILD CACHE** (e.g., Docker build volumes like `gnina-build-cache`).

If you encounter build errors like CUDA version mismatches or stale object files:
- Touch/modify the specific source file to force recompilation of that file only
- Use `touch filename.cu` or similar to update the timestamp
- Do NOT use `docker volume rm` or clear cache directories
- Ask the user before taking any destructive action on build artifacts


# K8s Pod Reservations

When reserving K8s pods:
- **ALWAYS use the default 5 hours** unless the user explicitly specifies a different duration
- Do NOT override the default duration with shorter times


# Running Commands on K8s Pods

When running commands via `kubectl exec`, **ALWAYS redirect full output to a file on the remote pod**, then retrieve the output separately. This prevents kubectl connection timeouts from losing output for long-running commands.

**Pattern:**
```bash
# Run command, redirect stdout+stderr to file on pod
~/.claude/run_cli.py run "kubectl --context=west1 exec POD_NAME -n development -- bash -c 'COMMAND > /tmp/output.log 2>&1'"

# Then retrieve output
~/.claude/run_cli.py run "kubectl --context=west1 exec POD_NAME -n development -- tail -100 /tmp/output.log"
```

**Example:**
```bash
# Run benchmark, save output on pod
~/.claude/run_cli.py run "kubectl --context=west1 exec my-pod -n development -- bash -c 'cd /workspace && PYTHONPATH=. .venv/bin/python benchmark.py --flag > /tmp/benchmark.log 2>&1'"

# Check if done (look for final output lines)
~/.claude/run_cli.py run "kubectl --context=west1 exec my-pod -n development -- tail -30 /tmp/benchmark.log"
```

- Use unique log filenames if running multiple commands (e.g., `/tmp/bench_grid.log`, `/tmp/bench_mlp.log`)
- For very long commands, check progress with `tail` periodically


# GNINA Build Workflow

## Building GNINA

**ALWAYS use the local Docker build script:**
```bash
./build_docker_incremental.sh
```

- This builds gnina using the Docker build cache volume
- For incremental builds after code changes, just run the script again
- If ninja says "no work to do", touch the modified source file first: `touch gninasrc/main/main.cpp`
- **NEVER build directly on K8s pods** - always build locally with Docker

## Deploying to K8s Pods

The locally-built binary may have library mismatches with K8s pod environments (e.g., libtorch symbols). To deploy properly:

1. Use a K8s pod with the **same image** as the build: `us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:latest`
2. The pod should mount or have access to the same library versions
3. If you get symbol lookup errors like `undefined symbol: _ZNK3c105Error4whatEv`, it's a libtorch version mismatch


# Debugging Flyte Workflows

## Finding Compounds1D Parquet Paths

When a Flyte task fails with compound/docking errors, find the input parquet file to inspect the error column:

1. **Get the workflow metadata:**
   ```bash
   kubectl get flyteworkflows.flyte.lyft.com <execution-id> -n development -o json
   ```

2. **Find the data directory** from the workflow status:
   - Look at `status.nodeStatus.<node-id>.outputDir` for the failing node
   - Example: `gs://opta-gcp-rezotx-uc-us-central1/metadata/propeller/<project>-<domain>-<execution>/...`

3. **Download an inputs.pb file** from the nested node structure:
   ```bash
   gsutil cp gs://.../fau5gk1i/data/0/dn0/0/dn*/inputs.pb /tmp/inputs.pb
   ```

4. **Search for parquet URIs** in the protobuf (they're in `gs://rezo-flyte/scratch/serializable/`):
   ```python
   import re
   with open('/tmp/inputs.pb', 'rb') as f:
       content = f.read().decode('utf-8', errors='ignore')
   parquet_uris = re.findall(r'gs://rezo-flyte/scratch/[^\s]+\.parquet', content)
   ```

5. **Download and inspect the parquet:**
   ```bash
   gsutil cp gs://rezo-flyte/scratch/serializable/<hash>.parquet /tmp/compounds.parquet
   ```
   ```python
   import pandas as pd
   df = pd.read_parquet('/tmp/compounds.parquet')
   print(df['error'].unique())  # Check error messages
   ```


# Pull Requests

When creating or updating a PR, always show the Graphite URL to the user (from `gt submit` output).

# Graphite

**NEVER run `gt restack`**. Restacking can cause conflicts and unintended changes across the branch stack. If `gt submit` fails because it wants a restack, stop and ask the user how to proceed.
