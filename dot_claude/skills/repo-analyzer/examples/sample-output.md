# Sample Analysis Output

Example output from analyzing a hypothetical full-stack application repository.

---

# Codebase Analysis: acme-webapp

## Summary
A full-stack web application for e-commerce, built with a Python backend (FastAPI) and React frontend. The codebase follows a monorepo structure with clear separation between frontend, backend, and shared components.

## Language Breakdown
| Language   | Files | Code Lines | Comments | Blank |
|------------|-------|------------|----------|-------|
| TypeScript | 127   | 18,450     | 1,230    | 2,100 |
| Python     | 89    | 12,340     | 2,100    | 1,500 |
| SQL        | 23    | 1,890      | 450      | 280   |
| YAML       | 18    | 620        | 85       | 90    |
| JSON       | 12    | 340        | 0        | 20    |
| Dockerfile | 3     | 85         | 25       | 15    |

**Total**: 272 files, 33,725 lines of code

## Module Structure
```
acme-webapp/
├── frontend/                   # React frontend application
│   ├── src/
│   │   ├── components/        # Reusable UI components
│   │   ├── pages/             # Route pages
│   │   ├── hooks/             # Custom React hooks
│   │   ├── services/          # API client services
│   │   ├── store/             # Redux state management
│   │   └── utils/             # Utility functions
│   ├── tests/                 # Frontend tests
│   └── package.json           # subproject
├── backend/                    # FastAPI backend
│   ├── app/
│   │   ├── api/               # API endpoints
│   │   ├── core/              # Core configuration
│   │   ├── models/            # SQLAlchemy models
│   │   ├── schemas/           # Pydantic schemas
│   │   ├── services/          # Business logic
│   │   └── utils/             # Utility functions
│   ├── tests/                 # Backend tests
│   ├── alembic/               # Database migrations
│   └── pyproject.toml         # subproject
├── shared/                     # Shared types and constants
│   └── types/
├── infra/                      # Infrastructure code
│   ├── terraform/             # Cloud infrastructure
│   ├── k8s/                   # Kubernetes manifests
│   └── docker/                # Docker configurations
├── scripts/                    # Development scripts
├── docs/                       # Documentation
└── .github/                    # GitHub Actions workflows
```

## Key Concepts

### Domain
E-commerce platform handling:
- Product catalog management
- User authentication and authorization
- Shopping cart and checkout flow
- Order processing and fulfillment
- Payment integration (Stripe)

### Architecture
- **Style**: Monorepo with separate frontend/backend
- **Backend**: FastAPI with SQLAlchemy ORM, PostgreSQL database
- **Frontend**: React 18 with TypeScript, Redux Toolkit for state
- **API**: REST with OpenAPI documentation
- **Auth**: JWT-based authentication with refresh tokens

### Entry Points
- `backend/app/main.py` - FastAPI application entry
- `frontend/src/main.tsx` - React application entry
- `scripts/dev.sh` - Development environment startup
- CLI: `backend/app/cli.py` - Admin CLI commands

### Key Modules

**Backend**
- `api/` - REST endpoints organized by resource (users, products, orders)
- `services/` - Business logic layer (OrderService, PaymentService)
- `models/` - Database models with relationships
- `core/` - Configuration, security, and middleware

**Frontend**
- `components/` - Design system components (Button, Card, Modal)
- `pages/` - Route-based page components
- `store/` - Redux slices for cart, user, products
- `services/` - API client with TypeScript types

## Dependencies & Layers

```
┌─────────────────────────────────────┐
│           Frontend (React)          │
│  pages → components → hooks/store   │
└──────────────────┬──────────────────┘
                   │ REST API
┌──────────────────▼──────────────────┐
│           Backend (FastAPI)         │
│   api → services → models → DB      │
└──────────────────┬──────────────────┘
                   │
┌──────────────────▼──────────────────┐
│         PostgreSQL Database         │
└─────────────────────────────────────┘
```

**Dependency Direction**:
- Frontend depends on backend API contracts
- API layer depends on services
- Services depend on models and external integrations
- No circular dependencies detected

## Notable Patterns

1. **Repository Pattern**: Data access abstracted through repository classes
2. **Dependency Injection**: FastAPI's `Depends()` for service injection
3. **Feature Slices**: Redux organized by feature (cart, user, products)
4. **Shared Types**: TypeScript types generated from OpenAPI spec
5. **Database Migrations**: Alembic with auto-generated migrations
6. **Environment Config**: Pydantic Settings for configuration management
7. **Testing Strategy**:
   - Unit tests with pytest/vitest
   - Integration tests for API endpoints
   - E2E tests with Playwright
