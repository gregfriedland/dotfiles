# Language-Specific Patterns Reference

Detailed reference for identifying module boundaries and project structure across different programming languages.

## Python

### Boundary Markers
- `__init__.py` - Package marker (can be empty)
- `pyproject.toml` - Modern project configuration (PEP 517/518)
- `setup.py` - Legacy project setup
- `setup.cfg` - Declarative setup configuration
- `requirements.txt` - Direct dependencies

### Standard Structure
```
project/
├── src/
│   └── package_name/
│       ├── __init__.py
│       └── module.py
├── tests/
├── pyproject.toml
└── README.md
```

### Monorepo Patterns
- Look for multiple `pyproject.toml` files in subdirectories
- Namespace packages (`pkg.subpkg` structure)
- `src/` layout vs flat layout

## JavaScript / TypeScript

### Boundary Markers
- `package.json` - Project configuration and dependencies
- `tsconfig.json` - TypeScript configuration
- `jsconfig.json` - JavaScript configuration

### Standard Structure
```
project/
├── src/
│   ├── components/
│   ├── utils/
│   └── index.ts
├── tests/
├── package.json
└── tsconfig.json
```

### Monorepo Patterns
- `workspaces` field in root `package.json`
- `packages/` or `apps/` directories
- Lerna (`lerna.json`)
- Turborepo (`turbo.json`)
- Nx (`nx.json`, `workspace.json`)

### Framework Detection
Check `package.json` dependencies for:
- `next` - Next.js
- `react` - React
- `vue` - Vue.js
- `@angular/core` - Angular
- `express` - Express.js
- `fastify` - Fastify
- `nestjs` - NestJS

## Go

### Boundary Markers
- `go.mod` - Module definition
- `go.sum` - Dependency checksums

### Standard Structure
```
project/
├── cmd/
│   └── myapp/
│       └── main.go
├── internal/
│   └── pkg/
├── pkg/
│   └── library/
├── go.mod
└── go.sum
```

### Conventions
- `cmd/` - Main applications (each subdirectory is a separate binary)
- `internal/` - Private packages (cannot be imported externally)
- `pkg/` - Public library code
- Single `main.go` at root for simple projects

## Rust

### Boundary Markers
- `Cargo.toml` - Package manifest
- `Cargo.lock` - Dependency lock file

### Standard Structure
```
project/
├── src/
│   ├── lib.rs      # Library crate root
│   ├── main.rs     # Binary crate root
│   └── bin/        # Additional binaries
├── tests/          # Integration tests
├── benches/        # Benchmarks
├── examples/       # Example code
└── Cargo.toml
```

### Workspace Pattern
```toml
# Root Cargo.toml
[workspace]
members = ["crate1", "crate2"]
```

## Java

### Boundary Markers
- `pom.xml` - Maven configuration
- `build.gradle` - Gradle build (Groovy)
- `build.gradle.kts` - Gradle build (Kotlin)
- `settings.gradle` - Multi-project settings

### Standard Structure (Maven)
```
project/
├── src/
│   ├── main/
│   │   ├── java/
│   │   └── resources/
│   └── test/
│       ├── java/
│       └── resources/
└── pom.xml
```

### Multi-Module Pattern
```
project/
├── module-a/
│   └── pom.xml
├── module-b/
│   └── pom.xml
└── pom.xml (parent)
```

## Ruby

### Boundary Markers
- `Gemfile` - Dependencies
- `*.gemspec` - Gem specification
- `Rakefile` - Task definitions

### Standard Structure
```
project/
├── lib/
│   └── my_gem/
│       └── version.rb
├── spec/
├── Gemfile
└── my_gem.gemspec
```

### Rails Structure
```
rails_app/
├── app/
│   ├── controllers/
│   ├── models/
│   └── views/
├── config/
├── db/
└── Gemfile
```

## C# / .NET

### Boundary Markers
- `*.csproj` - Project file
- `*.sln` - Solution file
- `Directory.Build.props` - Shared build properties

### Standard Structure
```
solution/
├── src/
│   └── MyProject/
│       └── MyProject.csproj
├── tests/
│   └── MyProject.Tests/
│       └── MyProject.Tests.csproj
└── MySolution.sln
```

## PHP

### Boundary Markers
- `composer.json` - Dependencies and autoloading
- `composer.lock` - Dependency lock

### Standard Structure
```
project/
├── src/
├── tests/
├── vendor/
├── composer.json
└── composer.lock
```

### Framework Detection
Check `composer.json` for:
- `laravel/framework` - Laravel
- `symfony/framework-bundle` - Symfony

## Swift

### Boundary Markers
- `Package.swift` - Swift Package Manager manifest
- `*.xcodeproj` - Xcode project
- `*.xcworkspace` - Xcode workspace

### Standard Structure (SPM)
```
package/
├── Sources/
│   └── MyLibrary/
├── Tests/
│   └── MyLibraryTests/
└── Package.swift
```

## Common Patterns Across Languages

### Test Directories
- `tests/`, `test/`, `spec/`
- `__tests__/` (JavaScript convention)
- `*_test.go` files (Go convention)
- `test_*.py`, `*_test.py` (Python convention)

### Documentation
- `docs/`, `doc/`, `documentation/`
- `README.md`, `README.rst`, `README.txt`
- `CHANGELOG.md`, `HISTORY.md`

### Configuration
- `config/`, `configs/`, `conf/`
- `.env`, `.env.example`
- `settings/`

### CI/CD
- `.github/workflows/` - GitHub Actions
- `.gitlab-ci.yml` - GitLab CI
- `.circleci/` - CircleCI
- `Jenkinsfile` - Jenkins
- `.travis.yml` - Travis CI

### Infrastructure
- `terraform/`, `infra/`
- `k8s/`, `kubernetes/`
- `docker/`, `Dockerfile`
- `ansible/`, `playbooks/`
