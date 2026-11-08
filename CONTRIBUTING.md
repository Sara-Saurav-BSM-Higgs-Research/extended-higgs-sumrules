# Contributing to Extended Higgs Sector Sum Rules

This repository is a collaborative scientific software project developed by Saurav Bania and Sara Chopra.

Our objective is to create transparent, verifiable and reusable computational methods for analytical and phenomenological studies of extended Higgs sectors.

## 1. Development Philosophy

All contributions should prioritize scientific correctness, reproducibility, clarity and maintainability.

We use object-oriented programming where it provides meaningful abstractions, and pure functions for symbolic transformations and calculations where appropriate.

The main Python package is located under `src/higgs_research/`. Exploratory notebooks should call reusable library functions rather than duplicate implementations.

## 2. Development Environment

The initial supported Python version is 3.12.

To prepare a local development environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

The Python package is installed in editable mode, allowing changes to the source code to be reflected without reinstalling after every source edit.

When dependencies change in `pyproject.toml`, the environment must be updated through the installation command.

## 3. Code Quality

Python code should follow the project's Ruff configuration and use descriptive names, type annotations and appropriate docstrings.

Before submitting changes, run:

```bash
ruff check .
ruff format --check .
python -m pytest -v
```

Code formatting can be applied with:

```bash
ruff format .
```

New nontrivial functionality should include automated tests.

## 4. Scientific Verification

New analytical calculations should document:

- Relevant Lagrangians and physical assumptions
- Field definitions and normalization conventions
- Metric signature and hypercharge convention
- Coupling and mixing-angle definitions
- Symbolic derivation methodology
- References to relevant equations and publications
- Independent validation where feasible

We distinguish between symbolic identities, independently verified physical derivations and agreement with previously published expressions.

Discrepancies with source theses or other references should be documented and investigated.

## 5. Git Workflow

The `main` branch is the integrated research branch.

Substantial changes should be developed on dedicated branches, such as:

- `feature/2hdm-potential`
- `feature/gm-couplings`
- `feature/scattering-amplitudes`
- `feature/higgstools-integration`

Pull requests should summarize the physics objective, implementation, assumptions, verification and potential limitations.

Both collaborators should review substantial scientific changes whenever practical.

## 6. Jupyter Notebooks

Notebooks are used for interactive derivations, numerical exploration, visualization and explanation.

Reusable functionality belongs in the Python package.

Notebooks should execute from a fresh kernel without depending on hidden variables or previous interactive sessions.

## 7. Research Data

All external datasets must have documented provenance, versions and applicable usage conditions.

Large generated datasets should not be committed to Git by default.

Credentials, access tokens and other secrets must never be committed.

## 8. Authorship and Citation

Scientific and software contributions should be attributed accurately.

Git commits and pull requests document development history, while scientific publication authorship is determined by substantive contributions and agreement among collaborators.

The project will provide machine-readable citation metadata through `CITATION.cff`.
