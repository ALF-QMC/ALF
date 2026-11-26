# Template for `HOOKS.md` on `master_ai`

*Copy the part below the line to `HOOKS.md` in the root of `master_ai` (once), then add one section per feature.*

**Format, which `tools/hook_lint.py` relies on:**
- each hook starts with a heading `### <path of the baseline file>`;
- one line of explanation follows;
- then the exact added lines, in a fenced code block.

`hook_lint.py` compares these blocks with `git diff master` and reports anything that differs.

The example below shows a feature `EXAMPLE` with a dispatch hook in `Wrapgr_mod.F90` and its build hooks.

---

# Hooks of the AI features in baseline files

*`master_ai` only. Written with AI assistance; not part of the ALF baseline.*

A hook is the only permitted edit to a file that exists on `master` (`AI_development/AI_DEVELOPMENT_RULES.md`, R2).
The diff of every baseline file against `master` consists of the blocks below and nothing else.

## Feature EXAMPLE (one-line description)

Code: `Prog/Example/`. Tests: `testsuite/Example.tests/`. Switch: `configure.sh ... EXAMPLE`.

### Prog/Wrapgr_mod.F90
Module `Wrapgr_mod`: makes the dispatch function visible.
```
#ifdef ALF_EXAMPLE
  use Example_mod, only: Example_wrapgrup
#endif
```

### Prog/Wrapgr_mod.F90
`WRAPGRUP`, first executable statement: the feature handles the slice, or returns `.false.` when it is off.
```
#ifdef ALF_EXAMPLE
    if (Example_wrapgrup(GR,NTAU,PHASE)) return
#endif
```

### Prog/Makefile
Adds the feature's sources when the switch is set.
```
ifneq (,$(findstring -DALF_EXAMPLE,$(ALF_FLAGS_PROG)))
SRCS += Example/Example_mod.F90
endif
```

### configure.sh
Reset of the switch (the script is sourced, so a value from an earlier call must not survive).
```
EXAMPLE_ENABLED=""
```

### configure.sh
Keyword.
```
    EXAMPLE)
      EXAMPLE_ENABLED="1"
    ;;
```

### configure.sh
Flag.
```
if [ "${EXAMPLE_ENABLED}" = "1" ]; then
  ALF_FLAGS_PROG="${ALF_FLAGS_PROG} -DALF_EXAMPLE"
fi
```

### Documentation/doc.tex
User page of the feature.
```
\input{example}
```
