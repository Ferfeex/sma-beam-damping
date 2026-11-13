# Changes from the original thesis code

The original code of the thesis was a collection of MATLAB scripts: the same solver existed as six near-identical
copies (one per case study), with hard-coded constants, large result files and several older attempts next to the
final method. This file lists what was kept, moved, removed and changed, and how the equivalence with the original
was checked. Bugs and oddities of the original are in `KNOWN_ISSUES.md`; the speed-up is in `PERF_NOTES.md`; the map of
the original code is in `docs/CODE_MAP_original_analysis.md`.

## Result in one paragraph

The final method (`Loss/FEM1D.m` + `myfun.m` + `fun_ode_1D.m` of the original) is one solver (`src/sma_beam_solve.m`)
driven by parameter files (`cases/case_*.m`, each returns a struct) through one entry point, `run_case('base')`.
This code reproduces the original: the short runs of every case equal the unmodified original scripts to round-off
(`tests/run_tests.m`), and the base case gives the mid-span ring-down peaks of thesis Fig. 4.3.

## Layout

| | original | here |
|---|---|---|
| copies of the base solver | 6 `FEM1D.m` + 6 `myfun.m` variants | 1 `sma_beam_solve.m`, 2 output functions (`myfun`, `myfun_SE`) |
| duplicates | 20 byte-identical copies of helper scripts | one copy of each |
| large data | `.mat` files of up to 168 MB, 24 `.fig` files | not versioned (`.gitignore`); only the small inputs `data/p.mat`, `data/p_SE.mat` and two reference files for the tests |
| encoding / line endings | 14 files Windows-1252, all CRLF | UTF-8, LF |

### Kept (main path)
* **`src/`** - the final method: `sma_beam_solve.m`, `myfun.m` (Newton on stress and martensite fractions + stress-reversal
  restart, the physics), `myfun_SE.m`, `fun_ode_1D.m`, `build_beam.m` (the FEM block that was pasted into every script),
  `linear_beam_solve.m` (linear reference and linear pulses), `composite_matrices.m`, `initial_velocity.m`, `meshare.m`.
* **`cases/`** - `case_base`, `case_hle_linear` (the linear reference), `case_composite`, `case_pulse` (the forced case),
  `case_prestrain_se`, `case_linear_pulse_incremental`, `case_linear_pulse_rectangular`, and `default_config.m` listing every
  parameter (the ones a menu fills are marked `[MENU]`: load type, force, frequency `f_ref`, duration `Tf`, initial velocity
  `init.scale`, number of nodes `geom.N`).
* **`brinson/`** - the material-level scripts that run in Octave: `Tension_compression_parametric.m` (the loss-factor fit),
  `PrestrainSE.m`, `Brinson.m`, `Brinson_validation.m`, `MyModel.m`, `Phasediagram.m`, `SuperElasticity.m`,
  `TC_parametri_paper.m`, `Red_ampl_range.m`. The Italian comments are unchanged.
* **`data/`** - the original loss-factor polynomials `p.mat` (default) and `p_SE.mat`.
* `compat/octave/ode23t.m` (stand-in, Octave only), `src/smooth_mean2d.m` (own moving average replacing a third-party
  smoothing function that came without a licence).

### Moved to `archive/` (not maintained, see `archive/README.md`)
Older attempts (`ODE_NEWT/`, `Mod/`, root DAE files, `FEM2D/`, case variants `Furbo_Pre` and `xiT0nullo`, broken Brinson
helpers) and `archive/original_case_scripts/`, the unmodified originals of what became cases, so that every change can
be diffed. Only `.m` sources were archived.

### Removed
* The byte-identical duplicates (one copy of each kept), all large data and results, the `.fig` files.
* `Loss/KFF_Computation.m` (the solver builds `KFF` itself; the hand-saved `KFF.mat` does not match it, see KNOWN_ISSUES #7)
  and other by-products of earlier ports (`octave_*.mat`, logs, `octave-workspace`).

## Changes in the code

| Change | Why | Numerics |
|---|---|---|
| 6 copies of the FEM script -> `sma_beam_solve` + case structs | one place to fix; menu-ready | identical (tested) |
| `myfun` takes a struct `prm` instead of 36 positional arguments (original names unpacked at the top) | readability | identical |
| hard-coded `750` and `X(31:60)` -> `numel(xisp)`, `X(2*N-1:end)` | N-independent code | identical for N = 16 |
| hard-coded reference frequency `/(2*pi*10)` (also 14.5, 18, 19, 21) in two files -> `cfg.damping.f_ref` | it had to be edited by hand in two places | identical |
| optional stiffness/damping of the GFRP skins (`KFF_ADD`, `CFF_ADD`; zeros for the pure SMA beam) | merges the composite variant of `myfun`; adding exact zeros changes no bits | identical |
| optional force pulse inside `myfun` / `fun_ode_1D` (`FFF`) | merges the transient variant | identical |
| globals (`p1`, `m1`, `tt`, `yy`, ...) are reset at the start of every solve | the original relied on `clear all` at the top of the script | fixes re-runs in the same session |
| restart-loop safety cap (`max_restarts`); progress printed every 0.05 s of simulated time | avoid endless loops and 4800-line logs | none |
| `Tension_compression_parametric.m`: `eps_upp = 0.001:0.001:0.025`, `eps_L = 0.016/-0.0146`, saves `results/p_regen.mat`, never overwrites `data/p.mat` | the stored version could not produce `p.mat` | affects only the regenerated `p` |
| `meshare.m`: section thickness and length are arguments; the title shows the real time | the 30 x 20 mm section and `t = 0.8` were hard-coded | plotting only |
| Newton loop condition `&&` -> `||` (`prm.newton_maxit = 100` caps the loop) | see KNOWN_ISSUES #1 | peaks change by <= 0.006 mm, curve by <= 0.8 % |
| `smooth2a.m` -> `smooth_mean2d.m` (conv2 moving average) | third-party code without licence | identical to 4e-16 (plot only) |
| mid-span read-out at exactly x = L/2 (cubic Hermite from nodes 8 and 9) instead of node 9 (x = 533 mm) | N = 16 has no node at L/2; see KNOWN_ISSUES #9 | `res.midspan` is 0.56 % larger for mode 1; the node value is kept as `res.midspan_node` |
| speed-up (vectorised assembly, analytic Jacobian, no per-step history, ...) | runtime | rounding level, see `PERF_NOTES.md` |

**Not changed:** all physics and numerics of the main method, variable names, Italian comments.

## The `&&` finding (details in KNOWN_ISSUES.md #1)

In the original `myfun.m` the Newton loop is `while ~isempty(indp) && ~isempty(indm)`: it only runs while there are points in
both a tension and a compression transformation zone and stops as soon as either set has converged; `||` looks like the
intent and is used here. Full 0.8 s base case, `&&` versus `||`: peaks 7.087 / 5.494 / 4.155 / 3.076 / 2.239 / 1.611 versus
7.087 / 5.490 / 4.149 / 3.075 / 2.241 / 1.612 mm (node-9 read-out); largest difference of the whole mid-span curve 0.058 mm
(0.8 %); `||` needs 4883 instead of 4320 steps. The effect is negligible for the thesis conclusions. To reproduce the original
results exactly, change `||` back to `&&` inside the parentheses in `src/myfun.m` / `src/myfun_SE.m`.

## Verification

| Check | Result |
|---|---|
| Short runs against the unmodified original scripts, run in a child Octave (`tests/run_tests.m`; the copy of the original is patched only in the Newton `while` line): base Tf = 0.03, composite Tf = 0.03, pulse Tf = 0.025 (covers the ramp and the 0.02 s release), pre-strain SE Tf = 0.02 | equal to round-off: relative difference < 1e-5 required, about 1e-8 observed |
| Composite beam: matrices recomputed by `composite_matrices.m` against the matrices saved by hand in 2020 | agree to 1e-9 relative (the 2020 file was written by MATLAB's `eig`) |
| Same case twice in one Octave session | identical (globals are reset) |
| Full base case (0.8 s) | peaks at x = 500 mm 7.13 / 5.52 / 4.17 / 3.09 / 2.25 / 1.62 / 1.16 / 0.82 mm; thesis Fig. 4.3 shows about 5.5, 4.15, 3.1, 2.2, 1.6 mm after the first |
| Full composite beam (0.8 s) | peaks 7.09 / 6.28 / 5.53 / 4.84 / 4.21 / 3.65 mm |
| Regenerated loss-factor fit against the original `data/p.mat` | agrees within 0.001-0.003 for strain >= 0.3 %, lower at small strain (0.119 vs 0.137 at 0.06 %) |

Runtimes in Octave (1 core): see the table in `README.md` and `PERF_NOTES.md`.

**Not verified:** full-length runs of `pulse` and `prestrain_se` against the original (hours; only shortened runs), anything in
MATLAB (`ode23t` in MATLAB is a different integrator than the Octave stand-in), N other than 16, the `archive/` code.
