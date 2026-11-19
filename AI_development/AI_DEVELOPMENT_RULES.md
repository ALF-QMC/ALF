# Rules for AI-assisted development in ALF

*Maintained by the ALF e.V. Changes to this file are made by humans, through a normal review on `master`.*

These rules apply to everybody who writes code for ALF with the help of an AI assistant (Claude, Copilot, ChatGPT, …), for any feature.

## The four basic rules (ALF e.V., 2026-11-19)

1. Never commit AI-written code into `master`.
2. Never change the baseline, i.e. the code path of `master`.
3. The baseline stays human-readable and keeps ALF's structure. The rules below define what that means in practice.
4. AI-written code lives on the branch `master_ai`; see `BRANCHING.md`. AI code reaches `master` only by a human rewrite.

This folder (`AI_development/`) is human-maintained documentation and tooling. It lives on `master` and is not part of the code path.

## Definitions

| Term | Meaning |
|---|---|
| **baseline** | Every file that exists on `master`, and the code path that a build of `master` executes. |
| **AI feature** | A self-contained extension written with AI assistance. Each feature has one short name `<FEAT>`, e.g. `GPU` or `LOWRANK`. |
| **AI switch** | The build switch that turns a feature on: a `configure.sh` keyword `<FEAT>` that adds `-DALF_<FEAT>`. Without it, the feature is not compiled. |
| **hook** | The only permitted edit to a baseline file: a small, listed block that hands control to an AI feature. The format is in R2. |

---

## R1. Where AI code lives

1. **All AI logic goes in new files.** Each feature gets its own subdirectories:
   - `Prog/<Feature>/` for QMC-core code;
   - `Libraries/Modules/<Feature>/` for model-independent libraries;
   - `testsuite/<Feature>.tests/` for its tests.
2. **Files follow ALF naming.**
   - `<Name>_mod.F90` for modules, `<Name>_smod.F90` for submodules.
   - One module per file.
   - Public names carry the module prefix, as in the baseline (`Op_*`, `Hop_mod_*`, `Predefined_*`).
3. **Use a Fortran submodule when an AI feature needs private data of a baseline module,** instead of making that data public. The parent module then gets one `interface` hook. ALF already uses submodules for its Hamiltonians.
4. **Runtime parameters of an AI feature go in their own namelist** (e.g. `&VAR_<FEAT>`), read and broadcast by the feature's own module.
   - The precedent is the tempering namelist `/VAR_TEMP/`, read by `read_and_broadcast_TEMPERING_var` in `QMC_runtime_var_mod.F90`.
   - `VAR_QMC` and `QMC_runtime_var_mod.F90` stay untouched.
   - If the namelist is missing from `parameters`, the feature is off.
5. **A new model goes in `Prog/Hamiltonians/`** and is registered through a file in `Prog/Hamiltonians.list.d/`. This drop-in mechanism needs no hook.

## R2. Hooks: the only edits allowed in baseline files

1. **A Fortran hook is at most three lines.** The middle line is either a `use` of the feature's module or a call into the feature, typically a "dispatch and return":
   ```fortran
   #ifdef ALF_<FEAT>
         if (<feat>_routine(...)) return
   #endif
   ```
   - No `#elif` or `#else`.
   - The baseline code below the hook stays exactly as on `master`.
   - The feature returns `.false.` when it is switched off at run time, so the baseline code then runs.
2. **Prefer hooks inside kernel routines over hooks in `main.F90`.** This is ALF's own pattern for variants:
   - STAB1/2/3/LOG branch inside `WRAPUR`, `WRAPUL`, `CGR` and `UDV_State%decompose`;
   - the callers never see the variant;
   - a hook in a kernel also covers all its other callers.
3. **For hopping, add a variant as a class.** `OpT_postprocess` in `Hop_mod.F90` already chooses between `RealExpOpT`, `CmplxExpOpT` and `OpT_time_dependent`. A further `ContainerElementBase` subclass selected there costs one hook.
4. **Every hook is listed in `HOOKS.md` on `master_ai`**, as the exact block of added lines, under the file it belongs to; see `HOOKS_TEMPLATE.md`. CI checks this list with `tools/hook_lint.py`. Each feature has a hook budget, agreed in review; keep it small.
5. **Build and documentation files** get longer blocks. These also count as hooks and are listed:
   - `configure.sh`: one keyword arm per feature, following the `HDF5` pattern, plus one flag block. Also one reset line for the switch variable, because the script is sourced.
   - `Prog/Makefile` and `Libraries/Modules/Makefile`: one block per feature that extends `SRCS` when the switch is set.
   - `Documentation/doc.tex`: one `\input`.
6. **If a feature needs a hook point that the baseline does not offer** (for example, an inline block in `main.F90` that a three-line hook cannot reach), ask the maintainers for an *extension point*: a human refactors `master` without changing its output, for example by moving the block into a named routine. Do not work around it with longer hooks.

## R3. The baseline path is untouched without the switch

1. A build without any AI switch must produce **bit-identical output** to `master` with the same compiler, flags and seeds. `tools/baseline_guard.py` checks this.
2. The same must hold for a build **with** the switch when the feature is switched off at run time.
3. Without the switch, nothing of the feature is compiled: no modules, no namelist reads, no environment-variable reads, no extra output in `info`.
4. **Corrections or speedups of baseline code** are not made inside baseline files on `master_ai`. Either they live in the feature, used only with the switch, or they are proposed for a human rewrite on `master`.

## R4. Style of AI files

AI files follow the baseline wherever the baseline is consistent, and are stricter where it is not.

| Topic | Rule |
|---|---|
| Licence | ALF GPL header with the "Section 7 additional terms" |
| Doxygen | Each module and public routine has a `!>` block with `@brief`, and `@param` for every argument |
| Declarations | `Implicit none` in every module and procedure; `Kind=Kind(0.d0)` and `cmplx(…, kind(0.d0))` |
| Language | Fortran 2008. It must compile with GNU, Intel, IntelLLVM and PGI/NVHPC in `Devel` mode (`-std=f2008 -fcheck=all`) |
| Visibility | `private` by default, explicit `public` |
| Errors | Message to `error_unit`, then `Terminate_on_error(ERROR_<X>, __FILE__, __LINE__)` (`runtime_error_mod.F90`). No `stop`, `error stop` or private exit wrappers |
| Warnings | `error_unit` with a `Warning:` prefix; repeats suppressed |
| info file | Open with `newunit`, rank 0 only, one line per setting the user can choose |
| Format | Spaces, consistent indentation within a file, comments in English |

**Forbidden in AI files:**
- Debug, self-test, fault-injection or timing code inside production modules. It belongs in `testsuite/<Feature>.tests/`, or in one separate `<Feature>_debug_mod.F90` behind its own macro.
- Environment variables that change numerics or behaviour. All runtime choices go through the feature's namelist. Variables read by vendor runtimes rather than by ALF belong in job scripts and on the feature's documentation page.
- References to tickets, chat sessions, agents, private notes or decision numbers. Use plain reasons instead.
- Commented-out code, `#line` directives, and padding `#ifdef` blocks.
- Code copied into several places. Use one routine.

## R5. Numerics and documentation

1. **Every place where an AI feature computes differently from the baseline** is listed on the feature's user page, with the reason and the measured effect. Examples: a different pivoting, a different summation order, an extra re-orthogonalisation.
2. **Every AI feature has one user page** in `Documentation/` (LaTeX, included from `doc.tex`). It covers:
   - the build switch;
   - the supported and refused cases;
   - the namelist parameters;
   - the differences in numerics;
   - how the feature was verified.
3. **Each merge into `master_ai` gets a changelog entry** in `CHANGELOG_AI.md` on `master_ai`, in the format of `CHANGELOG.md` (`### YYYY-MM-DD title`, Author, list). `CHANGELOG.md` itself is a baseline file.
4. **`AI_FEATURE_STATUS.md` on `master_ai` is updated in the same commit** as every feature that is ported, verified, refused or dropped. It is a feature × platform table; the format is in `FEATURE_STATUS_TEMPLATE.md`.

## R6. Tests

1. **Each feature has tests** in `testsuite/<Feature>.tests/`. Fortran unit tests are registered in CTest like those in `Prog.tests/`.
2. **Each feature is compared with the baseline:**
   - deterministically at kernel level (fixed fields or a fixed seed), with a documented tolerance;
   - at observable level, either by bin-by-bin comparison on the same seed (when the feature reproduces the Markov chain) or statistically (within 3σ, with χ² reported) against a build of the baseline with the same compiler.
3. **A test is shown to fail on a planted error** before it is used as a gate.
4. **Compare builds only with identical compiler flags.** A `Devel` build is not bit-identical to an optimised build.

## R7. CI on `master_ai`

| Job | Checks |
|---|---|
| `lint_stopcommands` | unchanged from `master` |
| `hook_lint` | `tools/hook_lint.py`: the diff against `master` in baseline files is exactly the blocks listed in `HOOKS.md` |
| `baseline_guard` | `tools/baseline_guard.py`: build `master` and `master_ai` without switches (GNU, noMPI, `OMP_NUM_THREADS=1`), run the serial scenarios of `testsuite/test_branch_parameters.json` with a fixed `seeds` file, and compare all output byte by byte, ignoring `info` and `info_git` |
| test matrix | the existing jobs, once without and once with each AI switch the runner can build |

`testsuite/compare_dirs.py` is **not** enough for the baseline guard: it compares analysed results with `np.allclose`, not bit by bit.

## R8. Provenance and the way back to `master`

1. Commits on `master_ai` that contain AI-written code name the assistant in a commit trailer, e.g. `Co-Authored-By: …`.
2. Each AI file says in its Doxygen header: `Written with AI assistance on master_ai; not part of the ALF baseline.`
3. **Code reaches `master` only by a human rewrite.**
   - A human writes the code on a normal `master` branch and may cite the `master_ai` commit as a reference.
   - `master_ai` is never merged into `master`, and nothing is cherry-picked from it.

## Working with an AI assistant: practical advice

- Give the assistant this file and `BRANCHING.md` at the start of a session. For Claude Code, refer to them from a `CLAUDE.md` in your working copy.
- Let it work in its own worktree on an `ai/feat/<feature>-<topic>` branch, never in a checkout of `master`.
- Before every commit, run `tools/hook_lint.py` and `tools/baseline_guard.py`; see `tools/README.md`.
- Treat the assistant's claims that something is "tested" or "equivalent" as claims until a check you can rerun shows it.

---

## Appendix: ALF's structure

The conventions that the rules above build on, as found on `master` in November 2026.

**Directories**
- `Prog/`: QMC core, one `*_mod.F90` per component, plus `main.F90`.
- `Prog/Hamiltonians/`: one submodule per model, plus a template.
- `Prog/Hamiltonians.list.d/`: drop-in model registration.
- `Libraries/Modules/`: model-independent utilities, built into `modules_90.a`.
- `Libraries/libqrref/`: vendored LAPACK.
- `Analysis/`: post-processing.
- `Documentation/`: the LaTeX manual (`doc.tex`).
- `testsuite/`: CTest unit tests, `test_vs_ed/`, `test_branch_*`.

**Build**
- `configure.sh` is one keyword `case`. Any unknown word becomes `MACHINE`, so a misspelled switch is silently taken as the machine name.
- The Makefiles have explicit `SRCS` lists; `gen_deps.py` derives the dependencies.
- Extra compiler flags come from `ALF_FLAGS_EXT`.

**How optional features plug in today**
- Through preprocessor flags in `main.F90`: MPI, TEMPERING, PARALLEL_PARAMS, HDF5.
- As flag-guarded calls into their own module:
  - tempering (`Exchange_Step`);
  - global moves (`Global_Updates`);
  - Langevin/HMC/MALA (the object `Langevin_HMC`).
- As kernel variants that branch inside the kernel (STAB1/2/3/LOG).
- As hopping variants, which are classes chosen in `OpT_postprocess`.

**Conventions**
- Licence header; Doxygen `!> @brief`.
- `Implicit none`; `Kind(0.d0)`.
- Errors via `Terminate_on_error`, with the codes in `runtime_error_mod.F90`.
- `stop` is forbidden by CI outside `Analysis/`.

**Inconsistent in the baseline itself, so not demanded of AI code beyond R4:**
- capitalisation, indentation width and line length;
- module name vs file name;
- default visibility, GPL headers and Doxygen coverage;
- info unit 50 vs `newunit`.
