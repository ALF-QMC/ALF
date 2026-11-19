# Branching for AI-assisted development (`master_ai`)

*Maintained by the ALF e.V. It completes rule 4 of `AI_DEVELOPMENT_RULES.md`.*

## Picture

```
master ─────●──────────●──────────●───────────→   humans only
             \          \          \              sync: merge master → master_ai, never the reverse
master_ai ────M──────────M──────────M─────────→   baseline + AI features (new files + listed hooks)
              \    ↑      \     ↑
ai/feat/<x>    ●─●─┘       ●──●─┘                  short-lived AI work, off master_ai, merged back after review + CI
ai/proj/<x>    ────────M────────────M─────────→    long-lived project lines, merge only from master_ai
tags           ai/prod/<date>-<run>                code state of a production run (immutable)
               ai/archive/<name>                   finished or abandoned AI branches
```

## Branches and tags

| Name | Role | Who merges | Rules |
|---|---|---|---|
| `master` | baseline | humans only | Receives no AI commits. |
| `master_ai` | AI line = `master` + AI features | maintainers, after review | Receives `master` by merge; it is never rebased or force-pushed. Receives `ai/feat/*` by merge once CI is green. |
| `ai/feat/<feature>-<topic>` | one piece of AI work | the author | Branches off `master_ai` and merges back. Then it gets an `ai/archive/…` tag and is deleted. |
| `ai/proj/<project>` | a project line not meant for general use, e.g. physics variants of one paper | the project owner | Merges from `master_ai` only; it is never merged into `master_ai`. Anything general that comes out of it goes back as a new `ai/feat/*` branch. |
| `ai/prod/<date>-<run>` (tag) | code state of a production run | whoever submits the run | Set at submission. Production is built from the tag, never from a moving branch. Keep a table "dataset → tag" with the data. |
| `ai/archive/<name>` (tag) | keeps the last commit of a deleted AI branch | the author | — |

The `ai/` prefix keeps the AI line visibly apart from `master` and from human feature branches in the shared repository.

## Flows

1. **`master` → `master_ai`.** Merge every time `master` changes, or at least weekly.
   - Conflicts can arise only at hook lines, because AI code lives in new files.
   - The maintainer doing the merge resolves them on `master_ai`.
2. **AI work.** Branch `ai/feat/x` off `master_ai` → work → CI green → review → merge into `master_ai` → tag `ai/archive/feat/x` → delete the branch.
3. **`master_ai` → `master`: only by human rewrite.**
   - A human opens a normal branch off `master` and writes the feature by hand. They may cite the `master_ai` commit as a reference.
   - AI commits are never merged or cherry-picked into `master`.
   - When the human version lands on `master`, the AI version is removed from `master_ai` in the next `master` → `master_ai` merge, so nothing exists twice.
4. **Extension points.** A feature may need a hook point that the baseline cannot offer, for example a routine boundary inside `main.F90`. A human then creates it on `master` as a pure refactoring with identical output. The AI feature hooks into it on `master_ai`.

## CI gates on `master_ai` (every push, and every merge into it)

1. `lint_stopcommands`: unchanged.
2. `hook_lint`: `AI_development/tools/hook_lint.py`.
3. `baseline_guard`: `AI_development/tools/baseline_guard.py`, requiring bit-identical output without switches.
4. The existing test matrix, plus the same jobs with each AI switch the runner can build. GPU features need a GPU runner.

A merge into `master_ai` needs all four gates green and one maintainer review.

## Open questions for the e.V.

1. Does `master_ai` live in the main ALF repository (GitHub and the GitLab mirror) or in a separate repository?
2. Which machine runs the GPU tests, and with what budget?
3. Is `master_ai` released or cited at all? Proposal: no releases, only `ai/prod/*` tags for reproducibility.
4. Who may merge into `master_ai`: the core developers, or a named group?
5. How are AI files marked: by the header line and the commit trailer (R8), or also by a directory prefix?
