# AI-assisted development in ALF

This folder holds the rules and tools for developing ALF with AI assistants. The ALF e.V. maintains it; changes go through a normal review on `master`. Nothing in here is part of the ALF code path.

| File | What it is |
|---|---|
| `AI_DEVELOPMENT_RULES.md` | **The rules.** Read them before you start. |
| `BRANCHING.md` | Branches and tags for AI work: `master_ai`, `ai/feat/*`, `ai/proj/*`, `ai/prod/*`, plus the CI gates. |
| `HOOKS_TEMPLATE.md` | Format of `HOOKS.md` on `master_ai`, the list of every edit to a baseline file. |
| `FEATURE_STATUS_TEMPLATE.md` | Format of `AI_FEATURE_STATUS.md` on `master_ai`, the feature × platform status table. |
| `tools/` | Checks that every commit on `master_ai` must pass; see `tools/README.md`. |

## Quick start for a new AI feature

1. Branch `ai/feat/<feature>-<topic>` off `master_ai`, in its own worktree.
2. Give your assistant `AI_DEVELOPMENT_RULES.md` and `BRANCHING.md` as context.
3. Put all new code in `Prog/<Feature>/` or `Libraries/Modules/<Feature>/`, behind the build switch `<FEAT>`.
4. Enter baseline files only through hooks of at most three lines, and list each one in `HOOKS.md`.
5. Before each commit, run `tools/hook_lint.py` and `tools/baseline_guard.py`.
6. Add tests, a user page, an entry in `CHANGELOG_AI.md` and a row in `AI_FEATURE_STATUS.md`.
7. Open a merge request into `master_ai`. You need green CI and one maintainer review.
