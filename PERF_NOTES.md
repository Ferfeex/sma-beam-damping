# PERF_NOTES — speed-up of the solver

All numbers: GNU Octave 10.3 on a Linux machine with 8 cores, one core per run (several runs were sometimes executing at the
same time, so wall times are indicative, +-~30 %). Reference = the code before the speed-up (ode23t -> `ode23s` stand-in,
numerical Jacobian).

## What changed (all switchable: `cfg.perf.<flag> = 0` gives the original code path; with all flags 0 the result is bit-identical to the code before the speed-up)
| flag | change | files |
|---|---|---|
| `vec` | the 16 `trapz(csi,trapz(-y,..))` calls per element (x15 elements, x2 for C) become one matrix product with precomputed weights; `Extr'*Kk*Extr` assembly becomes a sparse scatter-add; strain `B*Extr*u` precombined | `src/assemble_element_matrix.m`, `myfun.m`, `myfun_SE.m`, `sma_beam_solve.m` |
| `hist` (=0) | `tt`, `yy`, and per-step `epsyylon/xysp/xysm/xyT/s` were appended at every step (`[a x]`, O(n^2) copying) and never read back unless `post.keep_fields` | same |
| `fast_rhs` | constant `inv(M)` instead of `M\` in every right-hand-side call | `src/fun_ode_fast.m` |
| `jac` | analytic Jacobian `[-M^-1 C, -M^-1 K; I, 0]` handed to the solver (ode23s otherwise builds it with 2x60 extra RHS calls per step) | `src/jac_ode_fast.m` |
| `sets` | `setdiff/union/ismember/unique` on index lists replaced by logical masks (same index sets) | `myfun.m`, `myfun_SE.m` |

The defaults are all fast flags on (`cases/default_config.m`). Solver, tolerances, event/restart logic, Newton, physics: untouched.

## Results (same Tf, same solver ode23s, default tolerances of each case)
Base 0.8 s: 433 s -> 9 s (~45x); prestrain_se 0.15 s: 19.7 s -> 0.9 s (~22x); pulse 0.05 s: 458 s -> 9.4 s (~49x). Same number of
steps/restarts, peaks equal to round-off (1e-9..1e-12) as long as the restart sequence is the same.

Timings of the fast path: base 0.8 s: 9 s; composite: 1.5 s; prestrain_se 0.8 s: 41 s; pulse 0.8 s: 112 s.

The optional recording of the martensite fractions at ~100 times (`cfg.post.n_snapshots`, used by the martensite maps and the
fields animation) copies 3 x 750 numbers about 100 times and does not change the solution (bit-identical); the base case takes the
same time with and without it.

## Equivalence
* **Results match to rounding level.** Short runs: same number of steps, differences 1e-9..1e-6 relative. Full base run (0.8 s): the 8 peaks differ from the code before the speed-up by at most **0.15 %** (3.0746 vs 3.0722 mm at peak 4; one restart flipped, 147 vs 157 restarts). They are not bit-identical, and cannot be, because the sums are done in a different order.
* **Super-elastic case:** the first ~5 peaks agree to 1e-6; late peaks (after ~0.3 s) are sensitive to round-off and to tolerances (a 1e-5 relative change of AbsTol moves them by ~2 %; 3.6 % between the old and the new code at 0.4 s, peak 8). This is a property of the model plus solver (restarts at stress reversals on the step grid), not of the speed-up.
* `tests/run_tests.m` compares with the original scripts with a relative tolerance of 1e-5.
* With `plot_fields` (`cfg.post.keep_fields`) only the LAST step of the fields is recorded (all the plots use it); `cfg.post.keep_field_history = true` records every step (slow, O(n^2)).

## Caveats
* **Results are not converged in the tolerance.** Base case (RelTol 1e-3 default): peak 8 = 0.8124 mm; RelTol 1e-4: 0.8290; 1e-5: 0.8311; 1e-6: 0.8336.
  The ring-down is therefore ~2 % lower with the default tolerance than with a tight one (restarts are detected on the step grid).
  The super-elastic case (RelTol 1e-2 as in the thesis) is far more sensitive: with RelTol 1e-3 the peak sequence is different.
  Changing AbsTol by 1e-5 relative changes late SE peaks by ~2 %.
* Because of this, any change that perturbs round-off (even `vec`) can flip one restart and move late peaks by 0.1 % (base) to a few % (SE, late peaks, 0.4 s+). `tests/run_tests.m` therefore uses a relative tolerance instead of strict equality.
* Octave `ode15s` (SUNDIALS IDA) works on base for Tf=0.1 (4x faster than the fast ode23s) but throws "IDASolve failed" at t~0.027 (Tf=0.8) / 0.082 (Tf=0.3) and on prestrain_se. Not usable.
* `ode45`/`ode23` are incompatible with the OutputFcn interface (called with several time points) and are explicit (stiff problem) - not pursued.
* MATLAB was not available: the analytic Jacobian option and the vectorised code use only standard functions, but MATLAB timing/behaviour was NOT measured.
