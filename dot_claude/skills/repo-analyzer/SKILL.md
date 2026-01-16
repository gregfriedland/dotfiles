---
name: repo-analyzer
description: This skill should be used when the user asks to "analyze a repository", "summarize a codebase", "what languages are used", "show me the module structure", "understand this project", "give me a codebase overview", "analyze the code", or mentions analyzing code organization, project architecture, repository structure, or codebase analysis.
version: 1.0.0
---

# Repository Codebase Analyzer

Analyze the structure, languages, and architecture of a codebase. Focus on the code itself, not repository metadata like commits, contributors, or stars.

## Purpose

Provide comprehensive codebase analysis including:
- Language breakdown with lines of code
- Module and package structure
- High-level concepts and domains
- Architecture patterns and layers
- Entry points and interfaces

## When to Use

Activate this skill when the user wants to:
- Understand a new codebase
- Get an overview of project structure
- Know what languages are used and in what proportion
- Identify modules and their responsibilities
- Understand the architecture of a project

## Analysis Workflow

### Step 1: Language and LOC Analysis

Run the LOC analysis script to get language breakdown:

```bash
~/.claude/skills/repo-analyzer/scripts/analyze_loc.sh /path/to/repo
```

Parse the JSON output to create a language breakdown table showing:
- Language name
- Number of files
- Lines of code (excluding blanks/comments)
- Comment lines
- Blank lines

If `scc` is not installed, the script falls back to `cloc`. If neither is available, use file extension analysis with the Glob tool to estimate languages present.

### Step 2: Module Structure Detection

Run the structure detection script:

```bash
python3 ~/.claude/skills/repo-analyzer/scripts/detect_structure.py /path/to/repo
```

This script identifies:
- Package boundaries (via `__init__.py`, `package.json`, `go.mod`, etc.)
- Subproject markers (via `pyproject.toml`, `Cargo.toml`, etc.)
- Standard directory conventions (`src/`, `lib/`, `tests/`, etc.)
- Monorepo patterns

If the script is unavailable, manually identify structure by:
1. Listing top-level directories
2. Looking for boundary marker files
3. Reading configuration files that define project structure

### Step 3: Concept Extraction

Analyze discovered modules to identify:

**Core Domains**
- What business logic or functionality does the code implement?
- What problem domain does it address?

**Infrastructure**
- Utilities and helper functions
- Configuration management
- CI/CD and deployment code
- Database and storage layers

**External Interfaces**
- REST/GraphQL APIs
- CLI tools
- Web UIs or frontends
- Library exports

**Testing Strategy**
- Test organization (unit, integration, e2e)
- Test frameworks used
- Coverage approach

**Design Patterns**
- Architecture style (monolith, microservices, serverless)
- Code organization (MVC, clean architecture, domain-driven)
- Common patterns (factory, singleton, observer, etc.)

### Step 4: Architecture Overview

Determine:
- **Dependency flow**: Which modules depend on which others
- **Layer structure**: Presentation, business logic, data access
- **Entry points**: Main files, CLI scripts, API endpoints, exports
- **Configuration**: How the application is configured

## Output Format

Present analysis in this structured format:

```markdown
# Codebase Analysis: [repo-name]

## Summary
[1-2 sentence description of what this codebase does]

## Language Breakdown
| Language   | Files | Code Lines | Comments | Blank |
|------------|-------|------------|----------|-------|
| [lang]     | [n]   | [n]        | [n]      | [n]   |

**Total**: [n] files, [n] lines of code

## Module Structure
[Tree diagram showing directory structure with annotations]

## Key Concepts

### Domain
[What the code does, the problem it solves]

### Architecture
[Architecture style and patterns used]

### Entry Points
[Main files, CLI commands, API endpoints]

### Key Modules
[Brief description of each major module/package]

## Dependencies & Layers
[Description of how modules relate to each other]

## Notable Patterns
[Any interesting patterns, conventions, or approaches observed]
```

## Script Reference

### analyze_loc.sh
Location: `~/.claude/skills/repo-analyzer/scripts/analyze_loc.sh`

Runs `scc` (preferred) or `cloc` (fallback) to count lines of code. Outputs JSON format for parsing.

Usage:
```bash
./analyze_loc.sh /path/to/repo
```

### detect_structure.py
Location: `~/.claude/skills/repo-analyzer/scripts/detect_structure.py`

Walks the directory tree and identifies:
- Package boundaries per language
- Subproject markers
- Standard conventions

Usage:
```bash
python3 detect_structure.py /path/to/repo [--max-depth N] [--json]
```

## Fallback Approach

If scripts are unavailable, perform analysis manually:

1. **Language detection**: Use Glob to find files by extension, estimate proportions
2. **Structure detection**: List directories, look for config files
3. **Concept extraction**: Read key files (README, main modules) to understand purpose

## Additional Resources

- See `references/language-patterns.md` for language-specific boundary markers
- See `examples/sample-output.md` for a complete analysis example

## Best Practices

1. Start with LOC analysis to understand scale and primary languages
2. Identify the project type (library, application, monorepo) early
3. Focus on understanding the "why" not just the "what"
4. Note any unusual patterns or deviations from conventions
5. Keep the summary concise but comprehensive
