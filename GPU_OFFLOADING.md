# Offloading BLAS/LAPACK to Intel GPUs: changes and findings

Branch `intelgpu-blas`. Nothing is committed yet.

## Summary

- `configure.sh` has a new option `GPU`. With ifx it compiles ALF with OpenMP offloading and oneMKL's GPU variants, e.g. `source configure.sh IntelLLVM MPI GPU`.
- Every BLAS/LAPACK call in `Prog/` and `Libraries/Modules/` now goes through a wrapper `ALF_<ROUTINE>`. A wrapper sends the call to the GPU only if the call is at least as large as the threshold of that routine.
  - The threshold is set at compile time (`-DALF_GPU_THRESH_<ROUTINE>=N`) or at runtime (environment variable `ALF_GPU_THRESH_<ROUTINE>=N`).
  - The default is `-1`, meaning never offload, so behaviour only changes once you set thresholds.
- A benchmark, `Benchmarks/gpu_linalg_bench.out`, measures the crossover size for each routine and checks GPU results against the CPU.
- **Correctness:**
  - CPU builds produce output bit-identical to `master`, and all 39 testsuite tests pass.
  - With the GPU, a full ALF run matches the CPU run to round-off (≤ 6e-13 relative).
  - All 25 GPU wrappers match the CPU to ≤ 5e-14.
- **Requirement:** ifx and oneMKL must come from the **same oneAPI release**. With ifx 2025.3 and MKL 2026.2, every offloaded call leaked about 3.4 MB of host memory. This exhausted the memory of the test machine.

## Usage

```sh
source configure.sh IntelLLVM MPI GPU     # requires ifx and MKLROOT
make
cd Benchmarks && make && ./gpu_linalg_bench.out 4096       # all routines up to N=4096
./gpu_linalg_bench.out 2048 ZGEMM ZGETRF                     # selected routines

export ALF_GPU_THRESH_ZGEMM=256           # offload ZGEMM calls with work >= 256^3
export ALF_GPU_OFFLOAD=0                  # disable offloading entirely
```

On rank 0, ALF prints the active device and thresholds at startup. `MKL_VERBOSE=2` shows which calls ran on the GPU: those lines end in `GPU0`.

### Threshold semantics

A threshold `N` is the size of an equivalent square problem. A call is offloaded if `N >= 0` and the call's work is at least:

| Routine group | Work of a call | Offloaded if work ≥ |
|---|---|---|
| GEMM, HEMM | `m·n·k` | `N³` |
| TRSM, TRMM | `m·n·(m or n)` | `N³` |
| GETRF, GEQRF | `m·n·min(m,n)` | `N³` |
| GETRS | `n²·nrhs` | `N³` |
| GETRI | `n³` | `N³` |
| UNGQR/ORGQR, UNMQR | `m·n·k` | `N³` |
| HEEV/SYEV, GESVD | `n³` or `m·n·min(m,n)` | `N³` |
| GEMV, GERU | `m·n` | `N²` |
| SCAL, COPY, SWAP, NRM2 | `n` | `N` |

This handles the rectangular shapes ALF uses (`Ndim × N_part`) through their work. The threshold values stay small enough for plain integer `-D` macros.

## Changes

### `configure.sh`
- New keyword `GPU` (help text, argument parsing, `GPU_ENABLED`).
- It refuses to run unless the compiler is ifx and `MKLROOT` is set.
- It removes `-static` (offloading needs the dynamically loaded GPU runtime) and `-qopenmp`.
- It adds `-fiopenmp -fopenmp-targets=spir64 -DALF_GPU -I$MKLROOT/include` and links with `-qmkl -fiopenmp -fopenmp-targets=spir64`. No explicit SYCL libraries are needed.
- New function `check_gpu`, which compiles and runs `check_gpu.f90` (one offloaded ZGEMM) in the temporary directory.
  - A compile or link failure is fatal.
  - A missing GPU only produces a warning, e.g. on a login node.

### Wrapper layer (`Libraries/Modules/`)
- **`alf_linalg_mod.F90`** holds the configuration:
  - thresholds with compile-time defaults;
  - `alf_linalg_init` (reads the environment variables, detects the device, prints the configuration);
  - `alf_linalg_set_threshold` and `alf_linalg_get_threshold`;
  - `alf_linalg_offload_active`.

  Initialisation is lazy, so Analysis and the tests need no changes. `Prog/main.F90` calls `alf_linalg_init` explicitly after `MPI_INIT` to print the configuration.
- **`alf_linalg_wrappers.F90`** contains the wrappers. Each has exactly the BLAS/LAPACK argument list.
  - **Offloadable (25):**

    | Kind | Routines |
    |---|---|
    | BLAS level 3 | ZGEMM, DGEMM, ZHEMM, ZTRSM, ZTRMM |
    | BLAS level 2 | ZGEMV, DGEMV, ZGERU |
    | BLAS level 1 | ZSCAL, ZCOPY, ZSWAP, DZNRM2 |
    | LAPACK | ZGETRF, DGETRF, ZGETRS, ZGETRI, DGETRI, ZGEQRF, DGEQRF, ZUNGQR, DORGQR, ZUNMQR, ZHEEV, DSYEV, ZGESVD |

  - **Pass-through only (9)**, because oneMKL has no OpenMP offload variant: ZGEQP3, ZGEEV, ZLACPY, DLACPY, ZLASET, ZLAPMR, ZLAPMT, ZLACRM, ZLARCM.
  - The offload path is `!$omp target data map(...)` + `!$omp dispatch` + the standard BLAS/LAPACK call. oneMKL's `declare variant` replaces the call with its GPU version.
  - Each argument is mapped exactly over the elements BLAS touches, `ld·(cols−1)+rows`. This is correct when callers pass sub-blocks like `A(5,7)` with `lda > m`, and never reads past the end of an array.
  - Offloaded LAPACK calls query and allocate their own device workspace, because GPU workspace sizes differ from CPU ones. Workspace queries (`lwork = -1`) always run on the CPU.
- **`Makefile`:** the two new files are added. Intel's `__genmod` interface files for the wrappers are deleted after compiling (see design notes).

### Call sites
- About 165 calls were renamed `CALL ZGEMM(` → `CALL ALF_ZGEMM(` and so on in:

  | Directory | Files |
  |---|---|
  | `Libraries/Modules` | `mymats_mod`, `Mat_subroutines_mod`, `entanglement_mod` |
  | `Prog` | `cgr1_mod`, `cgr2_2_mod`, `udv_state_mod`, `upgrade_mod`, `Operator_mod`, `UDV_WRAP_mod`, `wrapul_mod`, `Global_mod`, `WaveFunction_mod`, `QDRP_decompose_mod`, `Hamiltonians/LRC_mod` |

- `EXTERNAL ZGEMM/ZGETRF` were removed from `entanglement_mod`, and `DZNRM2` → `ALF_DZNRM2` in `udv_state_mod`.
- The only direct calls left are `ZGEQRF_REF` and `ZUNGQR_REF` into the vendored `libqrref`. `Libraries/libqrref`, `testsuite/` and the unbuilt `matTypes_mod.F90` are untouched.
- **Fixed in passing:** `QDRP_decompose_mod` declared `RWORK` for ZGEQP3 as COMPLEX instead of REAL. This was harmless, since the buffer was oversized.

### Other files
- `Benchmarks/` (new):
  - `gpu_linalg_bench.F90` and `Makefile`.
  - Sizes are N = 16·√2^k up to `NMAX`.
  - Times are medians of at least 3 repetitions and at least 0.2 s per measurement, and include the host↔device transfer of each call.
  - The GPU output is checked against the CPU output. Eigen- and singular values are compared, not vectors, and QR factors by the absolute values of R.
- Documentation:
  - `Documentation/running.tex`: `GPU` row in the configure table, plus a paragraph on thresholds and the benchmark.
  - `README.md`: short note on GPU offloading.
- `.gitignore`: `*.modmic` (written by ifx offload builds).

## Design notes

- **Wrappers are external procedures, not module procedures.** ALF often passes an element of an assumed-shape or pointer array as the start of a matrix or vector, e.g. `ALF_ZSCAL(N2, z, Mat(Op%P(I),1), N1)` with `Mat(:,:)`.
  - This works with implicit interfaces (as with BLAS today), but the standard forbids it with explicit interfaces. gfortran and ifx then fail with error #8284.
  - Rewriting such calls as array sections would break the strided calls.
  - So the wrappers keep the implicit-interface calling convention of BLAS. `ALF_DZNRM2` must be declared `REAL(Kind=Kind(0.d0))` at the call site.
  - For the same reason, the `__genmod` files that Intel's `-warn interfaces` generates for the wrappers are removed in the Makefile. Otherwise Intel Devel builds fail.
- **No round-robin device selection across MPI ranks.** On the test machine, device 1 is the integrated UHD 730. Device selection is left to `OMP_DEFAULT_DEVICE` / `ONEAPI_DEVICE_SELECTOR`.
- **Data is copied to and from the device for every call.** Keeping matrices resident on the device across calls, e.g. in `wrapul`/`wrapur`/`cgr`, is the natural next step. The wrapper layer is the place to add it.

## Verification

| Check | Result |
|---|---|
| GNU noMPI build vs `master`, standard `Start` run (6×6 Hubbard) | bit-identical outputs (only `info_git` differs) |
| `testsuite` (ctest, GNU noMPI) | 39/39 passed |
| Builds: GNU Devel (`-Werror`), ifx Devel, ifx noMPI+GPU, ifx MPI+GPU (`mpiifort -fc=ifx`) | all succeed |
| `configure.sh GNU noMPI GPU` | stops with a clear error |
| CI lint: no `stop`/`error stop` in `*.F90` | clean |
| Benchmark, all 25 GPU wrappers vs CPU | max. relative deviation 5e-14 |
| Short ALF run (2×2 sweeps, ifx 2025.3), all level-3/LAPACK thresholds 0 (1600 GPU calls) vs `ALF_GPU_OFFLOAD=0` | all observables agree to 2e-15 – 6e-13 relative |

A full ALF run with *every* routine offloaded was not completed. It hit the memory leak described below; it should be repeated with the matched toolchain. The level-1 and level-2 wrappers were verified in the benchmark only.

## Findings

### Toolchain: matching ifx and oneMKL versions are required
With **ifx 2025.3 and oneMKL 2026.2**, every `!$omp dispatch` into oneMKL leaked about **3.4 MB of host memory**. The cause was the combination of `libomptarget` from 2025.3 with `libsycl`/MKL from 2026.2.
- The leak happened regardless of routine and of the wrapper code: a bare `target data` region does not leak, and `target data` + `dispatch` does.
- An ALF run with full offloading was OOM-killed after about 11k calls. The OOM killer also terminated other processes on the machine.
- The same loop in C with icx 2026.2 did not leak.
- **With ifx 2026.2, 4000 dispatches cost 22 MB in total**, i.e. no leak. A newer runtime cannot simply be swapped in under an older binary: the ABI changed (`__tgt_get_interop_obj`).

### Offload overhead
The minimal cost of an offloaded call dropped from **3–14 ms** (ifx 2025.3) to **0.02–0.07 ms** (ifx 2026.2) for most routines. It remains 0.2–0.7 ms for ZGEQRF, DGEQRF and ZGESVD.

### Crossover sizes
Measured on an **Intel Arc B580** with an i5-14400 host (`OMP_NUM_THREADS=10`), ifx 2026.2 + oneMKL 2026.2, square problems N ≤ 4096, transfers included.

| Routine | Recommended threshold N | Speedup at N=4096 | Best speedup |
|---|---|---|---|
| DORGQR | 256 | 5.3× | 5.3× @ 4096 |
| ZGEMM | 256 | 2.5× | 2.6× @ 1448 |
| ZHEMM | 362 | 2.2× | 3.0× @ 2048 |
| DGEMM | 512 | 2.2× | 2.2× @ 4096 |
| ZTRMM | 512 | 3.0× | 3.0× @ 4096 |
| ZGETRS | 512 | 2.3× | 2.4× @ 2048 |
| ZTRSM | 724 | 2.1× | 2.3× @ 2048 |
| ZGETRI | 724 | 2.1× | 2.1× @ 4096 |
| ZUNMQR | 724 | 2.4× | 2.9× @ 2048 |
| DGETRI | 1024 | 3.5× | 3.5× @ 4096 |
| ZGETRF | 2048 | 1.6× | 1.6× @ 4096 |
| ZGEQRF | 2896 | 1.7× | 1.7× @ 4096 |
| DGETRF | 4096 | 1.2× | 1.5× @ 2048 |
| DGEQRF | 4096 | 1.8× | 1.8× @ 4096 |
| ZUNGQR | never | 0.93× | 0.94× |
| ZHEEV | never | 0.92× | 0.99× |
| DSYEV | never | 0.99× | 0.99× |
| ZGESVD | never | 0.95× | 0.98× |
| ZGEMV, DGEMV, ZGERU | never | ≈0.2× | ≈0.2× |
| ZSCAL, ZCOPY, ZSWAP, DZNRM2 | never | ≈0.01× | ≈0.05× |

Ready-to-use settings for this machine:

```sh
export ALF_GPU_THRESH_ZGEMM=256  ALF_GPU_THRESH_DGEMM=512  ALF_GPU_THRESH_ZHEMM=362
export ALF_GPU_THRESH_ZTRSM=724  ALF_GPU_THRESH_ZTRMM=512
export ALF_GPU_THRESH_ZGETRF=2048 ALF_GPU_THRESH_DGETRF=4096 ALF_GPU_THRESH_ZGETRS=512
export ALF_GPU_THRESH_ZGETRI=724 ALF_GPU_THRESH_DGETRI=1024
export ALF_GPU_THRESH_ZGEQRF=2896 ALF_GPU_THRESH_DGEQRF=4096
export ALF_GPU_THRESH_DORGQR=256 ALF_GPU_THRESH_ZUNMQR=724
```

What this means for ALF:
- The matrix-product routines (GEMM, TRSM/TRMM, UNMQR, GETRS/GETRI) pay off once `Ndim` is a few hundred.
- The QR factorisation used for stabilisation (ZGEQRF) only pays off for very large systems.
- Forming Q explicitly (ZUNGQR), eigen- and singular-value solvers, and all level-1/2 BLAS never pay off with per-call transfers.
- ZUNGQR, ZHEEV, DSYEV and ZGESVD give GPU results *exactly* equal to the CPU. That suggests MKL computes much of them on the host even in GPU mode, which is consistent with their speedup of about 1×.
- The thresholds depend on the GPU, the host and the number of MPI ranks sharing a GPU. Rerun the benchmark with the production configuration on the target machine.

## Open points

- Repeat the full ALF validation with every routine offloaded, using the matched ifx 2026.2 toolchain.
- Consider keeping matrices resident on the device across calls (see design notes). The QR and level-1/2 results suggest that per-call transfers limit what offloading can gain in the stabilisation path.
- The default thresholds in `alf_linalg_mod.F90` are all `-1`. Machine-specific defaults could be added to `configure.sh` machine entries via `ALF_FLAGS_EXT`.
- Benchmark rectangular shapes (`Ndim × N_part`) as used in ALF, in addition to square ones.
