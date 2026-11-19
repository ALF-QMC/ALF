# Template for `AI_FEATURE_STATUS.md` on `master_ai`

*Copy the part below the line to `AI_FEATURE_STATUS.md` in the root of `master_ai` (once). Then add rows as features arrive (rule R5.4).*

**What to fill in:**
- **Rows** are features or feature variants. A sub-row (`↳`) records a combination that a feature refuses, or a part of it that was left out on purpose.
- **Columns** are platforms. Use only the columns that apply to your features; a CPU-only feature needs only "CPU".
- **The `master_ai` column** gives the state on `master_ai`, the check that supports it, and the commit.

---

# Status of the AI features: feature × platform

*`master_ai` only. Updated in the same commit as every feature that is ported, verified, refused or dropped.*

## Legend

| Letter | Meaning |
|---|---|
| P | production use, with recorded checks |
| V | implemented and verified (same seed vs CPU, ED, or a scheme comparison) |
| I | implemented, verification not documented |
| O | opt-in, verified, not in production |
| R | refused at setup: the run stops with a message |
| NI | not implemented |
| D | dropped on purpose (reason in the row) |
| — | not applicable, or not there yet |

**Check levels:**
- **L0:** `hook_lint` plus the bit-identical `baseline_guard`.
- **L1:** same seed against a same-compiler CPU build, or a bin-by-bin comparison of variants.
- **L1b:** bitwise comparison with the switch on and off.
- **L2:** kernel test against LAPACK, with planted defects.
- **L3:** ED reference.
- **L4:** a production-size run with a restart test.

## Table

| Feature | CPU | NVIDIA | Intel | AMD | `master_ai` (state, check, commit) |
|---|---|---|---|---|---|
| Build without any AI switch = baseline | — | — | — | — | V: L0, `<commit>` |
| `<FEATURE>` | V | — | — | — | V: L1 (max rel. diff. …), compilers …, `<commit>` |
| ↳ `<FEATURE>` + `<option>` | R | — | — | — | R: the run stops; reason …, `<commit>` |

## History of this table

- `<date>` (`<commit>`): …
