# Research agent workflow

- Use Python 3.12 and SymPy; follow `pyproject.toml` and `docs/conventions/`.
- Preserve documented real-field, hypercharge, potential and normalization conventions.
- Use exact symbolic derivations and rational constants; state all assumptions.
- Never hard-code expected physics formulas into production calculations.
- Check mathematics independently, including signs, factors and special limits.
- Add fixed, independently checked numerical regression tests for new calculations.
- Keep changes small, type-annotated and within the human-approved task scope.
- Preserve existing scientific code, tests and unrelated user changes.
- Develop on GitHub feature branches; keep the current branch unless instructed otherwise.
- Human approval is required before merging; scientific acceptance belongs to humans.
- Never automate commits, pushes or PR merges.

## Orchestration

1. The main Codex session orchestrates; allow at most two concurrent subagents.
2. Assign approved implementation and tests to one `implementer`; no other code writer.
3. After implementation stops, have `physics_reviewer` and `code_reviewer` inspect the same results independently.
4. Reviewers must not modify files, spawn agents or coordinate answers with each other.
5. Give reviewers scope and conventions, but withhold the other reviewer's conclusions.
6. Reviewers return concrete discrepancies with file/line references and mathematical or test evidence.
7. The orchestrator runs full `python -m pytest`, `ruff check .` and `ruff format --check .` using existing dependencies.
8. Collect commands, outcomes, independent checks and unresolved issues in concise `notes/` reports.
9. Reviewers return findings to the orchestrator, who writes their reports; avoid lengthy terminal summaries.
10. Allow at most two implementer correction cycles, each followed by independent re-review and full checks.
11. Escalate unresolved issues after two cycles, or convention ambiguities earlier, to the human.
12. Passing checks or reviewer agreement does not authorize scientific acceptance or Git operations.
