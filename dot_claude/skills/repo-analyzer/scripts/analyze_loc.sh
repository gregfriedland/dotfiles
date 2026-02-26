#!/bin/bash
# analyze_loc.sh - Count lines of code by language
# Uses scc (preferred) or cloc (fallback) to analyze a repository

set -e

REPO_PATH="${1:-.}"

if [ ! -d "$REPO_PATH" ]; then
    echo "Error: Directory '$REPO_PATH' does not exist" >&2
    exit 1
fi

# Try scc first (faster, more features)
if command -v scc &> /dev/null; then
    scc --format json "$REPO_PATH"
    exit 0
fi

# Fall back to cloc
if command -v cloc &> /dev/null; then
    cloc --json --quiet "$REPO_PATH"
    exit 0
fi

# Neither tool available
echo "Error: Neither 'scc' nor 'cloc' is installed." >&2
echo "Install scc: brew install scc (or go install github.com/boyter/scc/v3@latest)" >&2
echo "Install cloc: brew install cloc" >&2
exit 1
