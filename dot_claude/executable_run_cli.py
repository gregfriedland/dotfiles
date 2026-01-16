#!uv run
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "bashlex",
# ]
# ///
"""Run CLI commands with output logging and whitelist validation.

This script runs CLI commands via subprocess, logging stdout/stderr to
timestamped files in /tmp/claude/output/. Commands must match a pattern
in ~/.claude/cmd_whitelist.txt to be executed. Uses bashlex to parse
compound commands (&&, ||, |, ;) and validates each subcommand.

Usage:
    run_cli.py run <command> [args...]
    run_cli.py run "command with args"
    run_cli.py run --dry-run <command>
    run_cli.py add <regex_pattern>
"""

import time
import argparse
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import bashlex

OUTPUT_DIR = Path("/tmp/claude/output")
# WHITELIST_FILE = Path.home() / ".claude" / "cmd_whitelist.txt"
PATTERNS = {
    "remove file (move to trash dir instead)": [re.compile("rm -[rf]"), re.compile("find .* -delete"), re.compile("find .* -exec rm -[rf] {} \\;")],
    "redirect (write to file instead)": [re.compile(">&"), re.compile(">>"), re.compile("<<"), re.compile("2>&1")],
    "subshell completion": [re.compile(r"\$\(")],
    "looping": [re.compile(r"for .* in .*; do .*; done")],
}
# SKIP_WHITELIST = True

# def load_whitelist() -> list[re.Pattern]:
#     """Load whitelist patterns from file.

#     Returns:
#         List of compiled regex patterns.
#     """
#     if not WHITELIST_FILE.exists():
#         return []

#     patterns = []
#     for line in WHITELIST_FILE.read_text().splitlines():
#         line = line.strip()
#         if line and not line.startswith("#"):
#             try:
#                 patterns.append(re.compile(line))
#             except re.error as e:
#                 print(
#                     f"Warning: Invalid regex pattern '{line}': {e}",
#                     file=sys.stderr,
#                 )
#     return patterns


# def has_file_redirect(command: str) -> bool:
#     """Check if command contains file redirects (>, >>, >&, etc.).

#     Allows stderr-to-stdout redirect (2>&1) but blocks file redirects.

#     Args:
#         command: The full command string.

#     Returns:
#         True if command contains file redirects, False otherwise.
#     """
#     found_redirect = False

#     def visit_node(node: bashlex.ast.node) -> None:
#         """Recursively visit AST nodes to find redirect nodes."""
#         nonlocal found_redirect
#         if node.kind == "redirect":
#             # Check if it's an output redirect (>, >>, >&, etc.)
#             if hasattr(node, "type") and node.type in (">", ">>", ">&"):
#                 # Allow 2>&1 (stderr to stdout) and >&2 (stdout to stderr)
#                 if node.type == ">&" and hasattr(node, "output"):
#                     if node.output in (1, 2):
#                         pass  # This is 2>&1 or >&2, allow it
#                     else:
#                         found_redirect = True
#                 else:
#                     found_redirect = True
#         if hasattr(node, "parts"):
#             for part in node.parts:
#                 visit_node(part)
#         if hasattr(node, "list"):
#             for item in node.list:
#                 visit_node(item)

#     try:
#         parts = bashlex.parse(command)
#         for part in parts:
#             visit_node(part)
#     except bashlex.errors.ParsingError:
#         # If parsing fails, do a simple regex check as fallback
#         # Allow 2>&1 and >&2 but block other redirects
#         cmd_cleaned = re.sub(r"2>&1", "", command)
#         cmd_cleaned = re.sub(r">&2", "", cmd_cleaned)
#         if re.search(r"(>|>>|>&)\s*\S", cmd_cleaned):
#             return True

#     return found_redirect


def extract_simple_commands(command: str) -> list[str]:
    """Extract all simple commands from a potentially compound command.

    Uses bashlex to parse the command and extract individual commands
    separated by &&, ||, |, ;, etc.

    Args:
        command: The full command string (may contain &&, ||, |, ;).

    Returns:
        List of individual command strings.
    """
    commands = []

    def visit_node(node: bashlex.ast.node) -> None:
        """Recursively visit AST nodes to find command nodes."""
        if node.kind == "command":
            # Extract the command text from the original string
            cmd_text = command[node.pos[0] : node.pos[1]]
            commands.append(cmd_text.strip())
        elif hasattr(node, "parts"):
            for part in node.parts:
                visit_node(part)
        elif hasattr(node, "list"):
            for item in node.list:
                visit_node(item)

    try:
        parts = bashlex.parse(command)
        for part in parts:
            visit_node(part)
    except bashlex.errors.ParsingError as e:
        # If parsing fails, treat the whole command as a single command
        print(
            f"Warning: Could not parse command with bashlex: {e}",
            file=sys.stderr,
        )
        commands.append(command)

    return commands if commands else [command]


def get_command_name(command: str) -> str:
    """Extract the command name (first token) from a command string.

    Args:
        command: A simple command string.

    Returns:
        The first token (command name) from the command.
    """
    parts = command.split()
    return parts[0] if parts else ""


def find_pattern_matches(command: str, patterns: list[re.Pattern]) -> set[str]:
    """Check if command matches any pattern.

    Args:
        command: A simple command string.
        patterns: List of compiled regex patterns.

    Returns:
        Set of patterns that match the command.
    """
    return set(pattern.pattern for pattern in patterns if pattern.search(command))


def validate_subcommands(
    full_command: str
    ) -> bool:
    """Validate that all subcommands in a compound command are allowed.

    Args:
        full_command: The full command string (may be compound).

    Returns:
        True if all subcommands are allowed, False otherwise.
    """
    subcommands = extract_simple_commands(full_command)

    for pattern_type, patterns in PATTERNS.items():
        pattern_matches = set()
        for subcmd in subcommands:
            pattern_matches.update(find_pattern_matches(subcmd, patterns))
        if len(pattern_matches) > 0:
            print(f"Error: The following disallowed patterns '{pattern_type}' were found: {', '.join(pattern_matches)}", file=sys.stderr)
            return False
    return True



def get_command_tokens(command: str) -> tuple[str, str]:
    """Extract first two tokens from command for filename.

    Args:
        command: The full command string.

    Returns:
        Tuple of (token1, token2) for use in filename.
    """
    parts = command.split()
    token1 = parts[0] if parts else "unknown"
    token2 = parts[1] if len(parts) > 1 else "cmd"

    # Clean tokens for use in filename
    token1 = re.sub(r"[^\w\-]", "_", os.path.basename(token1))
    token2 = re.sub(r"[^\w\-]", "_", token2)

    return token1[:20], token2[:20]


# def add_to_whitelist(pattern: str) -> int:
#     """Add a regex pattern to the whitelist file.

#     Args:
#         pattern: The regex pattern to add.

#     Returns:
#         0 on success, 1 on error.
#     """
#     # Validate the pattern
#     try:
#         re.compile(pattern)
#     except re.error as e:
#         print(f"Error: Invalid regex pattern '{pattern}': {e}", file=sys.stderr)
#         return 1

#     # Check if pattern already exists
#     if WHITELIST_FILE.exists():
#         existing = WHITELIST_FILE.read_text()
#         if pattern in existing.splitlines():
#             print(f"Pattern already exists in whitelist: {pattern}")
#             return 0

#     # Append to whitelist file
#     with open(WHITELIST_FILE, "a") as f:
#         f.write(f"{pattern}\n")

#     print(f"Added to whitelist: {pattern}")
#     return 0


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Run CLI commands with output logging and validation."
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # # 'add' subcommand
    # add_parser = subparsers.add_parser(
    #     "add",
    #     help="Add a regex pattern to the whitelist.",
    # )
    # add_parser.add_argument(
    #     "pattern",
    #     help="The regex pattern to add to the whitelist.",
    # )

    # 'run' subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Run a command.",
    )
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate the command without executing it.",
    )
    run_parser.add_argument(
        "command",
        nargs="+",
        help="The command to run (can be quoted or as separate args).",
    )

    return parser.parse_args()


def main() -> int:
    """Run the CLI command with logging."""
    args = parse_args()

    # Handle 'add' subcommand
    # if args.subcommand == "add":
    #     return add_to_whitelist(args.pattern)

    # Handle 'run' subcommand
    # Build command string
    if len(args.command) == 1:
        command = args.command[0]
    else:
        command = " ".join(args.command)

    # Unescape zsh's escaping of ! (history expansion character)
    command = command.replace("\\!", "!")

    # # Check for file redirects
    # if has_file_redirect(command):
    #     print(
    #         "Error: File redirects (>, >>, >&) are not allowed.\n"
    #         "Note: Output is already saved to a log file automatically by this "
    #         "wrapper, so redirects are unnecessary.",
    #         file=sys.stderr,
    #     )
    #     return 1

    # # Load and check whitelist
    # whitelist_patterns = load_whitelist()
    # if not whitelist_patterns:
    #     print(
    #         f"Error: No whitelist patterns found. Create {WHITELIST_FILE} "
    #         "with regex patterns (one per line).",
    #         file=sys.stderr,
    #     )
    #     return 1

    # Validate all subcommands in the command
    all_allowed = validate_subcommands(
        command
    )
    if not all_allowed:
        return 1
        # print(""
        #     "Error: The following subcommands are not in whitelist:\n",
        #     file=sys.stderr,
        # )
        # for cmd in disallowed:
        #     print(f"  - '{cmd}'", file=sys.stderr)

        # # Generate suggested pattern from first disallowed command
        # parts = disallowed[0].split()
        # token1 = parts[0] if parts else "cmd"
        # # Use second token only if it doesn't start with a dash (flag)
        # if len(parts) >= 2 and not parts[1].startswith("-"):
        #     suggested_pattern = f"^{token1} {parts[1]}.*$"
        # else:
        #     suggested_pattern = f"^{token1} .*$"

        # print(
        #     f"\nTo add to whitelist, run:\n"
        #     f"  uv run ~/.claude/run_cli.py add '{suggested_pattern}'",
        #     file=sys.stderr,
        # )
        # return 1

    # Dry-run mode: validate only, don't execute
    if args.dry_run:
        print(f"[DRY-RUN] Command validated: {command}")
        print("[DRY-RUN] Would execute but --dry-run specified. Exiting.")
        return 0

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Generate log filename
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    token1, token2 = get_command_tokens(command)
    log_file = OUTPUT_DIR / f"{timestamp}-{token1}-{token2}.log"

    # Run command
    print(f"Running: {command}")

    with open(log_file, "w") as f:
        t = time.perf_counter()
        process = subprocess.run(
            command,
            shell=True,
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True,
        )
        elapsed = time.perf_counter() - t

    # Count lines and print last 20 lines (truncated to 200 chars)
    with open(log_file, "r") as f:
        lines = f.readlines()
    line_count = len(lines)
    print(f"Log file: {log_file} (lines={line_count}, code={process.returncode}, elapsed={elapsed:.2f}s)")
    print("-" * 60)
    print("Last 20 lines:")
    for line in lines[-20:]:
        line = line.rstrip("\n")
        if len(line) > 200:
            line = line[:200] + "..."
        print(line)

    return process.returncode


if __name__ == "__main__":
    sys.exit(main())
