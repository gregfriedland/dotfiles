# Refactor Code

Refactor Python code following project style guidelines and design patterns.

## Usage
```
/refactor <file_or_directory> [instructions]
```

**Parameters:**
- `file_or_directory` - File path, directory, or glob pattern to refactor
- `instructions` (optional) - Specific refactoring goals or focus areas

## Examples

```bash
# Refactor a single file
/refactor rezo/atlas/library/scoring.py

# Refactor with specific goal
/refactor rezo/common/utils/gcs.py consolidate duplicate helpers

# Refactor a directory
/refactor rezo/atlas/library/enamine_bq/ extract shared logic into base class
```

## Instructions

When this command is invoked with arguments `$ARGUMENTS`:

1. Parse the arguments:
   - First argument is the file path, directory, or glob pattern
   - Remaining arguments are optional refactoring instructions

2. **Explore phase** - Before making any changes:
   - Read all target files fully
   - Search the codebase for duplicate classes/functions that overlap with the target code
   - Identify all callers/importers of public symbols being refactored
   - Identify all test files for the target code
   - Summarize what you found and what you plan to change. Ask the user to confirm before proceeding.

3. **Apply refactoring** following all rules below, then:
   - Run `ruff format` on every modified file
   - Run `ruff check --fix` on every modified file
   - Run `pyright` on every modified file
   - Fix any issues found
   - Run tests for affected code: `PYTHONPATH=$PWD pytest <test_files> -v`
   - Fix any test failures (fix the code, not the tests, unless the test is wrong)

4. **Report** what was changed and why.

---

## Refactoring Rules

### Structural Rules

- **Consolidate duplicates**: If a class or function already exists elsewhere with similar purpose, refactor to use the existing one. Extend it with additional fields/params if needed.
- **Move utility methods to domain classes**: Functions that operate primarily on a single class's data should be methods on that class.
- **Prefer method calls over standalone functions**: Replace `function(obj, args)` with `obj.method(args)` when the function logically belongs to the object's interface.
- **Required parameters**: Make previously optional parameters required if they are always needed.
- **Simplifying initialization**: Remove init parameters specific to certain implementations. Move them to factory methods or to where implementations are created.
- **Shared code location**: Move generally usable methods/functions to `rezo/common/utils/` or the appropriate shared utility directory rather than leaving them in domain-specific modules.
- **Colocate private methods with callers**: Private methods should be grouped in the same class as the public methods that call them. If a private helper is only used by one class, it belongs on that class, not as a standalone module function.
- **Remove unused code**: Delete functions, classes, imports, and constants that have no callers outside of tests. Search the codebase to confirm before removing.
- **Remove backwards-compatibility fallbacks**: Delete fallback code paths that exist only for backwards compatibility (e.g. old parameter names, deprecated enum values, try/except import shims for old library versions). This is an internal monorepo with no external consumers.
- **Resolve circular imports by restructuring**: NEVER use inline imports. Refactor the module structure instead, or ask the user how to proceed.

### Naming Rules

- **Descriptive method names**: Use names that clearly indicate purpose. E.g. `run_full_analysis` -> `run_analysis_all` to indicate it runs all types.
- **Descriptive class names**: Use specific names, not general concepts. E.g. `StatisticalTester` -> `QSetStatisticalTester` when specific to Q-set analysis.
- **File naming**: When renaming the main class in a file, rename the file too, including any test files.
- **Private naming**: Prefix internal functions/methods with `_`. E.g. `create_custom_library` -> `_create_custom_library`. Any module-level function or class method that is only called within the same file (not imported elsewhere) should be private, except test helpers.
- **Ignored variables**: Use `_` for variables that are unused.

### Class Design Rules

- **OOP**: Group related functions and state into classes. Use `__init__` to initialize state, `__str__` for a detailed yet fast string representation.
- **BaseModelNoExtra**: Classes with state should inherit from `rezo.common.utils.basemodel.BaseModelNoExtra` to enforce type safety (extra="forbid").
- **State + Logic together**: Group related configuration and methods in the same BaseModelNoExtra sub-class.
- **Abstract classes and factories**: Similar classes should inherit from a common (potentially abstract) base class. Use an enum with class names and add the factory method to the enum class. Generalize arguments so they apply to all classes.
- **Minimize public methods**: Make intermediate methods private with `_` prefix.
- **Functors**: Use functors when there is a clear high-level "action" method that calls private methods.
- **Member variables**: Initialize ALL member variables inside `__init__`, never in other methods. Convert member variables to local variables when they're only used in one method.
- **Enums**: Use `enum.Enum` for all int or string constants that are mutually exclusive. Never pass raw string value choices to functions.

### Size Rules

- **Function length**: Keep functions under 50 lines of code.
- **File length**: Keep files under 1000 lines of code.

### Style Rules

- **Imports**: Use absolute imports to the monorepo root. NEVER import inside functions/classes. Never use wildcard imports.
- **Line length**: Never exceed 80 characters (except long URLs with `# noqa: E501`).
- **Type hints**: All variables and returns for all methods. Use `list`, `dict`, `tuple` (not `List`, `Dict`, `Tuple`).
- **Callable types**: Use specific callable type hints like `Callable[[pd.DataFrame, pd.Series], np.ndarray]`, not generic `Callable`. Define a top-level type alias when used in multiple places.
- **Global constants**: Place all constants (public CAPS and private `_CAPS`) at the top of the file, after imports but before `logger` and any functions/classes. Public constants first, then private constants. No global variables.
- **Logging**: Create the logger after constants: `logger = logging.getLogger(__name__)`. No `print()` except in top-level scripts with `main()`.
- **Flyte v2 tasks**: Flyte v2 task functions must be decorated with `@rz_task` (from `rezo.common_v2.flyte.entity_v2`) and must always have a `_task` suffix (e.g. `_export_bq_to_gcs_task`, `cluster_enamine_library_task`).
- **Docstrings**: Google-style (pydocstyle) with `Args` and `Returns` sections. Avoid qualitative adjectives ("robust", "efficient", "comprehensive").
- **Exceptions**: Never swallow exceptions without at least debug-level output.
- **Init files**: Add blank `__init__.py` to new directories.
- **Type errors**: Add `# type: ignore` only for unexpected type errors.
- **GCS**: Use `from rezo.common.utils.gcs import GSPath` for GCS URLs.

### Design Pattern Rules

- **PydanticArgparse for CLI**: Auto-generate CLI from Pydantic models.
- **Nested config structure**: Nest configuration classes for hierarchical CLI args.
- **SerializableDataFrame**: Use for DataFrame fields inside Pydantic BaseModel classes (Flyte-serializable).
- **field_validator**: Use `@field_validator` in BaseModel classes for type conversions.
- **Consistent field naming**: When refactoring config classes, update all YAML/config files to match canonical class definitions.

### Config Pattern Example

```python
class ModelType(enum.Enum):
    RIDGE = "ridge"
    LASSO = "lasso"

class AnvilRunner(BaseModelNoExtra):
    data: SerializableDataFrame
    target_col: str
    model_type: ModelType

    def train(self) -> AnvilResults:
        ...

class AnvilCLIArgs(BaseModelNoExtra):
    anvil: AnvilRunner
    remote: bool = False
    flyte_project: FlyteProject
    flyte_domain: FlyteDomain

def main() -> None:
    parser = PydanticArgparse(AnvilCLIArgs)
    args = parser.parse_model()
    ...
```

### Code Organization Rules

- **Encapsulate library-specific logic**: When a module imports and uses a third-party library's API in scattered places, consolidate those usages into methods on the domain class that owns the configuration. Keep the orchestration module free of library internals.
- **Keep worker functions module-level**: Functions that need to be pickled for multiprocessing must stay at module level, not as methods. Place them in the same module as the class they serve.
- **Entry point at top of file**: Place the main orchestrator / public entry-point function as the first function in the file, above the private helpers it calls. Python resolves names at call time, so forward references are fine.
- **Merge identical constants**: When multiple constants have identical values, merge them into one with a descriptive name. Apply to config objects, node types, environment definitions, etc.
- **Inline single-use helpers**: If a private helper function is called from exactly one site and is short, inline it at the call site and remove the function to reduce indirection.
- **Propagate parameters through call chains**: When a function accepts an optional parameter (e.g. `columns`), ensure all code paths forward it. Don't silently ignore it for some file formats / branches.
- **Remove dead code near returns**: `del` statements right before a function `return` are useless since the frame is about to be destroyed. Keep `del` only when it frees large objects mid-function in memory-constrained contexts.
- **Simplify stale workarounds**: When you find defensive code added for a bug that's since been fixed (e.g. tuple-swap guards, type coercions), simplify to the straightforward version.

### Test Update Rules

- **Update test class names**: When renaming functions or moving methods, update test class names and docstrings to reflect new location/naming.
- **Remove trivial config tests**: Tests that only verify a config parameter is stored and returned (e.g. set `name="foo"`, assert `obj.name == "foo"`) don't test meaningful behavior. Keep tests that verify computations, transformations, or side effects.
- **Document future work**: Add TODO comments in docstrings for planned enhancements out of scope for the current PR.
