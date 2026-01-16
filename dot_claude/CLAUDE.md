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

You have access to a temporal knowledge graph via Graphiti MCP. Use it to:
- Remember facts about the user and project across sessions
- Track how information evolves over time
- Search for relevant context before making decisions

## MCP Tools Available

* mcp__graphiti__add_episode: Store new information (conversations, decisions, results)
* mcp__graphiti__search_facts: Find relevant relationships/edges before acting
* mcp__graphiti__search_nodes: Find entity summaries
* mcp__graphiti__get_episodes: Retrieve recent conversation history

## Protocol for using memory

* **ALWAYS CHECK THE MEMORY RIGHT AFTER THE USER SUBMITS A REQUEST:**
  1. Call `search_facts` with the user's query/topic. This returns edge info between nodes.
  2. Call `search_nodes` to get context on the nodes connected by relevant facts
  3. Incorporate relevant information into your context window
* **ALWAYS UPDATE THE MEMORY RIGHT BEFORE RESPONDING TO THE USER:**
  * Call `add_episode` with:
     - Key decisions made
     - Results/outcomes
     - User preferences learned
     - Important context for future
* If user says "Add X to memory": call add_episode

## What to Store (add_episode)

Parameters:
* name: identifier for the episode
* episode_body: content
* source_type | episode_type: enum: text (e.g. natural language description) or json (e.g. structure data from ml experiment)
* source_description: description of where the information came from (e.g. "ml experiment X" or "research paper Y")
* reference_time: current datetime
* group_id: always fixed to constant greg

## What NOT to Store

  - Raw data dumps or large outputs
  - Secrets, API keys, passwords
  - Temporary debugging info
  - Obvious/trivial information


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
