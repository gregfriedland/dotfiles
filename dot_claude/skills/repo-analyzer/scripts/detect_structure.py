#!/usr/bin/env python3
"""Detect module and package structure in a repository.

Identifies package boundaries, subproject markers, and standard directory
conventions across multiple programming languages.
"""

import argparse
import json
import os
from collections import defaultdict
from pathlib import Path

# Language-specific boundary markers
BOUNDARY_MARKERS = {
    "python": ["__init__.py", "pyproject.toml", "setup.py", "setup.cfg"],
    "javascript": ["package.json"],
    "typescript": ["package.json", "tsconfig.json"],
    "go": ["go.mod"],
    "rust": ["Cargo.toml"],
    "java": ["pom.xml", "build.gradle", "build.gradle.kts"],
    "ruby": ["Gemfile", "*.gemspec"],
    "php": ["composer.json"],
    "csharp": ["*.csproj", "*.sln"],
    "swift": ["Package.swift"],
    "elixir": ["mix.exs"],
    "haskell": ["*.cabal", "stack.yaml"],
}

# Standard directory conventions
STANDARD_DIRS = {
    "src": "Source code",
    "lib": "Library code",
    "pkg": "Packages",
    "cmd": "Command-line tools (Go)",
    "internal": "Internal packages (Go)",
    "app": "Application code",
    "tests": "Test files",
    "test": "Test files",
    "spec": "Test specifications",
    "docs": "Documentation",
    "doc": "Documentation",
    "examples": "Example code",
    "scripts": "Utility scripts",
    "bin": "Binaries/executables",
    "tools": "Development tools",
    "config": "Configuration",
    "configs": "Configuration",
    "assets": "Static assets",
    "public": "Public assets",
    "static": "Static files",
    "templates": "Template files",
    "migrations": "Database migrations",
    "fixtures": "Test fixtures",
    "mocks": "Mock implementations",
    "vendor": "Vendored dependencies",
    "node_modules": "Node.js dependencies",
    "dist": "Distribution/build output",
    "build": "Build output",
    "out": "Output directory",
    ".github": "GitHub configuration",
    ".circleci": "CircleCI configuration",
    "infra": "Infrastructure code",
    "deploy": "Deployment configuration",
    "k8s": "Kubernetes manifests",
    "terraform": "Terraform configuration",
    "ansible": "Ansible playbooks",
}

# Directories to skip
SKIP_DIRS = {
    ".git",
    ".svn",
    ".hg",
    "node_modules",
    "vendor",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    "venv",
    ".venv",
    "env",
    ".env",
    "dist",
    "build",
    "out",
    "target",
    ".next",
    ".nuxt",
    "coverage",
    ".coverage",
    "htmlcov",
    ".tox",
    "eggs",
    "*.egg-info",
}


def should_skip(path: Path) -> bool:
    """Check if a directory should be skipped."""
    name = path.name
    if name in SKIP_DIRS:
        return True
    if name.startswith(".") and name not in {".github", ".circleci"}:
        return True
    return False


def find_markers(directory: Path) -> dict[str, list[str]]:
    """Find boundary markers in a directory."""
    markers: dict[str, list[str]] = defaultdict(list)

    try:
        for item in directory.iterdir():
            if item.is_file():
                name = item.name
                for lang, lang_markers in BOUNDARY_MARKERS.items():
                    for marker in lang_markers:
                        if marker.startswith("*"):
                            if name.endswith(marker[1:]):
                                markers[lang].append(name)
                        elif name == marker:
                            markers[lang].append(name)
    except PermissionError:
        pass

    return dict(markers)


def get_dir_annotation(name: str) -> str | None:
    """Get annotation for a standard directory name."""
    return STANDARD_DIRS.get(name.lower())


def analyze_directory(
    root: Path,
    max_depth: int = 4,
    current_depth: int = 0,
    prefix: str = "",
) -> list[dict]:
    """Recursively analyze directory structure."""
    results = []

    if current_depth > max_depth:
        return results

    try:
        items = sorted(root.iterdir(), key=lambda x: (not x.is_dir(), x.name))
    except PermissionError:
        return results

    dirs = [item for item in items if item.is_dir() and not should_skip(item)]

    for i, item in enumerate(dirs):
        is_last = i == len(dirs) - 1

        markers = find_markers(item)
        annotation = get_dir_annotation(item.name)

        entry = {
            "name": item.name,
            "path": str(item.relative_to(root.parent)),
            "depth": current_depth,
            "is_last": is_last,
        }

        if markers:
            entry["markers"] = markers
        if annotation:
            entry["annotation"] = annotation

        # Check for subproject indicators
        if any(
            lang_markers
            for lang, lang_markers in markers.items()
            if any(
                m in lang_markers
                for m in [
                    "pyproject.toml",
                    "package.json",
                    "go.mod",
                    "Cargo.toml",
                    "pom.xml",
                ]
            )
        ):
            entry["is_subproject"] = True

        results.append(entry)

        # Recurse into subdirectories
        child_results = analyze_directory(
            item,
            max_depth=max_depth,
            current_depth=current_depth + 1,
            prefix=prefix + ("    " if is_last else "│   "),
        )
        results.extend(child_results)

    return results


def format_tree(results: list[dict], root_name: str) -> str:
    """Format results as a tree diagram."""
    lines = [f"{root_name}/"]

    depth_counters: dict[int, int] = defaultdict(int)
    depth_totals: dict[int, int] = defaultdict(int)

    # Count items at each depth
    for entry in results:
        depth_totals[entry["depth"]] += 1

    for entry in results:
        depth = entry["depth"]
        depth_counters[depth] += 1

        # Build prefix
        prefix_parts = []
        for d in range(depth):
            # Check if there are more items at this depth level
            if d < depth - 1:
                prefix_parts.append("│   ")
            else:
                if entry["is_last"]:
                    prefix_parts.append("└── ")
                else:
                    prefix_parts.append("├── ")

        prefix = "".join(prefix_parts)
        name = entry["name"]

        # Build annotation
        annotations = []
        if entry.get("annotation"):
            annotations.append(entry["annotation"])
        if entry.get("is_subproject"):
            annotations.append("subproject")
        if entry.get("markers"):
            langs = list(entry["markers"].keys())
            if langs and not entry.get("is_subproject"):
                annotations.append(f"{', '.join(langs)} package")

        annotation_str = f"  # {'; '.join(annotations)}" if annotations else ""

        lines.append(f"{prefix}{name}/{annotation_str}")

    return "\n".join(lines)


def detect_project_type(results: list[dict], root: Path) -> dict:
    """Detect overall project type and characteristics."""
    info: dict = {
        "type": "unknown",
        "languages": [],
        "is_monorepo": False,
        "subprojects": [],
    }

    # Collect all languages found
    languages = set()
    subprojects = []

    for entry in results:
        if entry.get("markers"):
            languages.update(entry["markers"].keys())
        if entry.get("is_subproject"):
            subprojects.append(entry["name"])

    info["languages"] = sorted(languages)
    info["subprojects"] = subprojects
    info["is_monorepo"] = len(subprojects) > 1

    # Determine project type
    root_markers = find_markers(root)
    if root_markers:
        primary_lang = list(root_markers.keys())[0]
        if "package.json" in root_markers.get("javascript", []):
            # Check for common frameworks
            pkg_path = root / "package.json"
            if pkg_path.exists():
                try:
                    with open(pkg_path) as f:
                        pkg = json.load(f)
                    deps = {
                        **pkg.get("dependencies", {}),
                        **pkg.get("devDependencies", {}),
                    }
                    if "next" in deps:
                        info["type"] = "Next.js application"
                    elif "react" in deps:
                        info["type"] = "React application"
                    elif "vue" in deps:
                        info["type"] = "Vue application"
                    elif "express" in deps:
                        info["type"] = "Express.js application"
                    else:
                        info["type"] = "Node.js project"
                except (json.JSONDecodeError, IOError):
                    info["type"] = "Node.js project"
        elif "pyproject.toml" in root_markers.get("python", []):
            info["type"] = "Python project"
        elif "go.mod" in root_markers.get("go", []):
            info["type"] = "Go project"
        elif "Cargo.toml" in root_markers.get("rust", []):
            info["type"] = "Rust project"
        elif root_markers.get("java"):
            info["type"] = "Java project"

    if info["is_monorepo"]:
        info["type"] = f"Monorepo ({info['type']})" if info["type"] != "unknown" else "Monorepo"

    return info


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Detect module and package structure in a repository"
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="Path to repository (default: current directory)",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=4,
        help="Maximum directory depth to analyze (default: 4)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output as JSON instead of tree diagram",
    )

    args = parser.parse_args()
    root = Path(args.path).resolve()

    if not root.is_dir():
        print(f"Error: '{args.path}' is not a directory", file=__import__("sys").stderr)
        __import__("sys").exit(1)

    results = analyze_directory(root, max_depth=args.max_depth)
    project_info = detect_project_type(results, root)

    if args.json:
        output = {
            "root": str(root),
            "project_info": project_info,
            "structure": results,
        }
        print(json.dumps(output, indent=2))
    else:
        print(f"Project Type: {project_info['type']}")
        if project_info["languages"]:
            print(f"Languages: {', '.join(project_info['languages'])}")
        if project_info["is_monorepo"]:
            print(f"Subprojects: {', '.join(project_info['subprojects'])}")
        print()
        print(format_tree(results, root.name))


if __name__ == "__main__":
    main()
