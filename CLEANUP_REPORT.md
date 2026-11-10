# CLEANUP REPORT

Cleanup of Federico Aghito's master's-thesis code (SMA beam vibration / damping) into a
repository-ready folder that runs in GNU Octave 10.3. Done on the night of 10 Nov 2026,
autonomously, on a copy; **nothing was pushed or published, there is no remote, and the
Google Drive original was never touched.**

Inputs (read-only): `/workspace/thesis-code/` (download of the Drive folder, 202 of 214 files) and
`/workspace/thesis-octave/` (earlier working Octave port, `OCTAVE_PORT_LOG.md`).

## 1. Result in one paragraph

The final method (`Loss/FEM1D.m` + `myfun.m` + `fun_ode_1D.m`, which existed as six near-identical
copies, one per case study) is now **one solver** (`src/sma_beam_solve.m`) driven by **parameter files**
(`cases/case_*.m`, each returns a struct) through **one entry point**, `run_case('base')`. On the full
base case the new code gives **exactly the same numbers as the unmodified original code** (max |difference|
= 0 over all 4320 time steps and 60 state variables, same Octave and integrator); for the other cases the
equivalence was verified on shortened runs (also exact, see section 5). The base run reproduces the
mid-span ring-down peaks **7.09 (t = 0), 5.49, 4.15, 3.08, 2.24, 1.61 mm** - peaks 2-6 are the
5.5 / 4.15 / 3.1 / 2.2 / 1.6 mm of thesis Fig. 4.3.

## 2. Before / after

| | before (thesis-code, local copy) | after (`sma-beam-damping`, git-tracked) |
|---|---|---|
| files | 202 (of 214 on the Drive) | 133 tracked files |
| size | 66 MB local (~996 MB on the Drive) | 0.9 MB (+ 1.8 MB `.git`, which holds the commit history and PNGs) |
| `.m` files | 114 (94 distinct, 20 byte-identical copies) | 113: 34 working code (incl. 8 case files, tests, tools), 65 archived older attempts, 14 unmodified originals kept for diffing (`archive/original_case_scripts`) |
| `.mat` files | 61 (up to 168 MB each) | 5 (all < 1 KB: `p.mat`, `p_SE.mat`, reference matrices, one archived `KFF.mat`+`p.mat`) |
| `.fig` | 24 | 0 (not versioned; thesis PDF has the plots; originals on the Drive) |
| copies of the base solver | 6 `FEM1D.m` + 6 `myfun.m` variants | 1 `sma_beam_solve.m`, 2 output functions (`myfun`, `myfun_SE`) |
| encoding / line endings | 14 files Windows-1252, all CRLF | UTF-8, LF |

## 3. What was kept, moved, removed, changed

### Kept (main path) - and why
* **`src/`** - the final method: `sma_beam_solve.m` (from `Loss/FEM1D.m`), `myfun.m` (Newton + stress-reversal
  restart, the physics), `myfun_SE.m`, `fun_ode_1D.m`, `build_beam.m` (the FEM block that was pasted into every script),
  `linear_beam_solve.m` (linear reference / linear pulses), `composite_matrices.m`, `initial_velocity.m`, `meshare.m`.
* **`cases/`** - `case_base`, `case_hle_linear` (the linear reference), `case_composite`, `case_pulse` (the
  `Loss/Transient` forced case), `case_prestrain_se` (`Loss/SE_pre-strain`), `case_linear_pulse_incremental` and
  `case_linear_pulse_rectangular` (`Transient_linear`), and `default_config.m` listing every parameter, with the
  fields a future menu would fill marked `[MENU]` (load type, force, frequency `f_ref`, duration `Tf`, initial
  velocity `init.scale`, number of nodes `geom.N`).
* **`brinson/`** - the Brinson material scripts that run in Octave: `Tension_compression_parametric.m` (the
  loss-factor fit), `PrestrainSE.m`, `Brinson.m`, `Brinson_validation.m` (**contains `%Tentativo 3 WINNER CAZZO`,
  still on line 76, byte for byte, file untouched except encoding/line endings**), `MyModel.m`, `Phasediagram.m`,
  `SuperElasticity.m`, `TC_parametri_paper.m`, `Red_ampl_range.m` (purpose unclear but it runs).
* **`data/`** - the author's original `p.mat` (default) and `p_SE.mat`.
* `compat/octave/ode23t.m` (stand-in, Octave only), `third_party/smooth2a.m` (with `THIRD_PARTY.md`).

### Moved to `archive/` (not maintained, see `archive/README.md`)
`ODE_NEWT/` (first Newton/fsolve approach), `Mod/` (DAE with `ode15s`), the root DAE files and scratch (`FEM_NL.m`,
`functxi.m`, `myfun2.m`, ...), `FEM2D/`, `Loss/Furbo_Pre` and `Loss/xiT0nullo` (case variants), the broken Brinson helpers, and
`archive/original_case_scripts/` = the unmodified originals of what was turned into cases (so each change can be diffed).
Only `.m` sources were archived (plus 2 tiny `.mat` of `xiT0nullo`).

### Removed
* **20 byte-identical duplicate scripts** (one copy of each was kept). Duplicate groups, with original paths:
- `Loss/Composite/KFF_Computation.m` = `Loss/KFF_Computation.m`
- `Loss/Composite/fun_ode_1D.m` = `Loss/SE_pre-strain/fun_ode_1D.m` = `Loss/fun_ode_1D.m` = `Transient_linear/Rectangular/fun_ode_1D.m` = `Transient_linear/fun_ode_1D.m`
- `Loss/Composite/meshare.m` = `Loss/Furbo_Pre/meshare.m` = `Loss/xiT0nullo/meshare.m`
- `Loss/Composite/smooth2a.m` = `Loss/Furbo_Pre/smooth2a.m` = `Loss/SE_pre-strain/smooth2a.m` = `Loss/Transient/smooth2a.m` = `Loss/smooth2a.m` = `Loss/xiT0nullo/smooth2a.m` = `ODE_NEWT/Variabili_1000N/smooth2a.m`
- `Loss/Furbo_Pre/fun_ode_1D.m` = `Loss/xiT0nullo/fun_ode_1D.m`
- `Loss/Furbo_Pre/myfun.m` = `Loss/xiT0nullo/myfun.m`
- `Mod/Condizioni_iniziali.m` = `Mod/DXI_algebric/Condizioni_iniziali_unforced.m`
- `Mod/DXI_algebric/fun_eps.m` = `Mod/fun_eps.m` = `ODE_NEWT/fun_eps.m`
- `Mod/fun_kff.m` = `ODE_NEWT/fun_kff.m`
- `Mod/meshare.m` = `meshare.m`
* **All large data and results** (`.mat` of 5-168 MB, ~1 GB on the Drive; the 12 files that exceeded the 10 MB download
  cap were never downloaded) and the 24 `.fig` files. `.gitignore` excludes `*.mat` (except the small inputs),
  `*.fig` and `results/`.
* `Loss/KFF_Computation.m` (the solver builds `KFF` itself; the hand-saved `KFF.mat` does not match it, see KNOWN_ISSUES #7),
  `Loss/octave_*.mat`, `full_run.log`, `octave-workspace` and other port by-products.

### Changed (with reasons)
| Change | Why | Numerics |
|---|---|---|
| UTF-8 + LF everywhere | Windows-1252/CRLF breaks diffs and GitHub display | none |
| 6 copies of the FEM script -> `sma_beam_solve` + case structs | one place to fix; menu-ready | identical (tested) |
| `myfun` takes a struct `prm` instead of 36 positional arguments (original names unpacked at the top) | readability | identical |
| hard-coded `750` and `X(31:60)` -> `numel(xisp)`, `X(2*N-1:end)` | N-independent code | identical for N = 16 |
| hard-coded reference frequency `/(2*pi*10)` (also 14.5, 18, 19, 21) in two files -> `cfg.damping.f_ref` | it had to be edited by hand in two places | identical |
| optional stiffness/damping of the GFRP skins (`KFF_ADD`, `CFF_ADD`; zeros for the pure SMA beam) | merges the Composite variant of `myfun`; adding exact zeros changes no bits | identical |
| optional force pulse inside `myfun` / `fun_ode_1D` (`FFF`) | merges the Transient variant | identical |
| globals (`p1`, `m1`, `tt`, `yy`, ...) are reset at the start of every solve | the original relied on `clear all` at the top of the script | fixes re-run-in-same-session |
| restart-loop safety cap (`max_restarts`) and progress printing every 0.05 s of simulated time (original printed every step) | avoid endless loops / 4800-line logs | none |
| `Tension_compression_parametric.m`: `eps_upp = 0.001:0.001:0.025`, `eps_L = 0.016/-0.0146`, saves `results/p_regen.mat`, never overwrites `data/p.mat` | the stored version could not produce `p.mat` | affects only the regenerated `p` |
| `meshare.m`: section thickness and length are arguments; title shows the real time | the 30 x 20 mm section and `t = 0.8` were hard-coded | plotting only |
| `Newton` loop condition: `&&` kept as default, `||` available as an opt-in flag | see KNOWN_ISSUES #1 | default identical |

**Not changed:** all physics and numerics of the main method, variable names, Italian comments.

## 4. The `&&` finding (details in KNOWN_ISSUES.md #1)

In `myfun.m` the Newton loop is `while ~isempty(indp) && ~isempty(indm)`: it only runs while there are
points in BOTH a tension and a compression transformation zone and stops when either set has converged; `||` looks
like the intent. Not silently changed. Full 0.8 s base case, `&&` versus `||`: peaks 7.087 / 5.494 / 4.155 / 3.076 /
2.239 / 1.611 versus 7.087 / 5.490 / 4.149 / 3.075 / 2.241 / 1.612 mm; max difference of the whole mid-span
curve 0.058 mm (0.8 %); `||` needs 4883 instead of 4320 steps (~13 % more). The effect is negligible for the thesis
conclusions. Federico decides whether to adopt it (`run_case('base','newton_or',true)`).

## 5. Verification

| Check | Result |
|---|---|
| **Full base case** (`run_case('base')`, Octave 10.3, 0.8 s) vs the stored reference run of the unmodified original | 4320 steps in both; **max abs difference of T and of all 60 state variables = 0** (bit-identical) |
| Base case peaks (mid-span, time axis from the first max) | 7.087 (t=0), 5.494 (0.1105 s), 4.155 (0.2203), 3.076 (0.3287), 2.239 (0.4371), 1.611 (0.5441), 1.147, 0.811 mm. Thesis Fig. 4.3: ~5.5, 4.15, 3.1, 2.2, 1.6 -> **match** |
| **Full composite** run (0.8 s, with the 2020 `Matrices.mat`) vs stored reference | 1153 steps both, **bit-identical** (difference 0) |
| Full composite run with the matrices recomputed by `composite_matrices.m` | different Rayleigh coefficients by ~4e-12 (the 2020 file was written by MATLAB's `eig`) -> the adaptive step sequence differs (1139 vs 1153 steps); mid-span curves differ by max 0.031 mm of 7.05 mm; peaks identical to 3 digits (7.054 / 6.243 / 5.497 / 4.810 / 4.187), 3.630 vs 3.629 mm at peak 6 |
| Linear reference (`hle_linear`) and linear incremental pulse vs the port's results | bit-identical (difference 0) |
| Short runs vs the **unmodified original scripts** run in a child Octave (`tests/run_tests.m`): base Tf=0.03 (581 steps), pulse Tf=0.025 (1912 steps, covers the ramp and the 0.02 s release), pre-strain SE Tf=0.02 (170 steps), composite Tf=0.03 | **all difference 0** (composite: 0 with the original `Matrices.mat`, 1.9e-9 relative with recomputed ones) |
| Same case twice in one Octave session | identical (globals are reset) |
| Regenerated loss-factor fit (`Tension_compression_parametric.m`) vs the original `data/p.mat` | regenerated `p` agrees with the original within 0.001-0.003 for strain >= 0.3 %, lower at small strain (0.119 vs 0.137 at 0.06 %) |
| `tests/run_tests('full')` (final state of the repo) | TESTS_FULL |

Runtimes in Octave (1 core): base ~370 s (333 s alone), composite ~53-80 s, `hle_linear` 2 s, linear pulses 10-17 s.
**Not verified:** full runs of `pulse` and `prestrain_se` (hours; only the shortened runs above), anything in MATLAB
(no MATLAB available; `ode23t` in MATLAB is a different integrator than the Octave stand-in), N other than 16, the
`archive/` code.

## 6. Open questions for Federico
1. **Licence**: none chosen (README says so). MIT or BSD-3-Clause? Is `smooth2a.m` (Greg Reeves 2009, no licence) acceptable in a public repo, or replace/drop it (only needed for the optional field plot)?
2. **`&&` vs `||`** in the Newton loop: keep the original behaviour (default) or switch? Effect ~0.8 % on the curve.
3. **p.mat**: keep the original hand-saved curve (default, matches the thesis best) or the regenerated one (lower at small strain)? Do you remember the amplitude set you used?
4. **Linear transient force term** was lost; is the inferred `M a = -K u - C v + F` right, and what were the exact loads of thesis Figs 3.21/3.23 (stored files: 50 N ramp, 37.5 N rectangular; the SMA `Transient` case has 200 N)?
5. Which archived variants (`Furbo_Pre`, `xiT0nullo`) correspond to which thesis figure (Fig. 3.20?) - should they become cases?
6. Do you want the original `.fig` files, thesis figures or the large result files in the repo (e.g. via Git LFS or a release), or only code?
7. Name/identity: the git commits use author `Federico Aghito <fed.aghito@gmail.com>` (his Gmail address, taken from the connected account). Change if you prefer another e-mail, especially before pushing (e-mail addresses in public commits are visible).
8. Repository name `sma-beam-damping` OK? README quotes the thesis title/author/year/supervisors from the PDF; remove names if you prefer.
9. `Red_ampl_range.m` (purpose unclear) and `Brinson.m` were kept in `brinson/` because they run; fine?
