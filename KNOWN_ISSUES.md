# Known issues of the original (2020) code

Found while reading, porting and cleaning the code. "Status" says what was done about each one.
Only #1 (the Newton loop condition) changes the numerical results, by < 1 %; everything else is documentation of problems in files that were archived or left as they were.

## 1. `&&` slip in the Newton loop (`src/myfun.m`, `src/myfun_SE.m`) - FIXED

```matlab
while ~isempty(indp) && ~isempty(indm)      % 2020 original
while (~isempty(indp) || ~isempty(indm)) && i < prm.newton_maxit      % now
```

`indp` = integration points in the tension transformation zone, `indm` = points in the compression zone.
With the original `&&` the Newton iteration **ran only while both sets were non-empty**, and **stopped as soon as
either set had converged**, even if the other had not. If only one set is active at a time step, the martensite
fractions of those points were not updated in that step. `||` ("iterate while any point still has to
converge") is the evident intent. **The 2020 code, and therefore the results produced for the thesis, used `&&`.**

**Status: corrected and made the default** (the `&&` option no longer exists; `prm.newton_maxit = 100` caps the
iterations so that a non-converging point cannot hang the loop). `tests/run_tests.m` compares the new code with
the original scripts after patching only this line in the copy of the original.

**Effect, full 0.8 s base case (same Octave, same integrator):**

| | original `&&` | `\|\|` (now) |
|---|---|---|
| solver steps | 4320 | 4883 |
| run time (indicative, machine load differed) | 370 s | 431 s |
| mid-span peaks [mm] | 7.087, 5.494, 4.155, 3.076, 2.239, 1.611, 1.147, 0.811 | 7.087, 5.490, 4.149, 3.075, 2.241, 1.612, 1.150, 0.814 |
| peak times [s] | 0, .1105, .2203, .3287, .4371, .5441, .6510, .7562 | 0, .1104, .2200, .3288, .4367, .5428, .6499, .7570 |
| max difference of the mid-span curve | | 0.058 mm (0.8 % of the 7.09 mm amplitude), rms 0.013 mm |

Peak heights differ by <= 0.006 mm; the thesis conclusions are unaffected, but `||` needs ~13 % more steps.
To reproduce the original numbers exactly, replace the `||` inside the parentheses by `&&` in `src/myfun.m` and `src/myfun_SE.m` (the results of the original
`&&` runs are quoted above).

## 2. `Condizioni_iniziali_unforced.m` uses an undefined variable (`XA`)
It defines `X0s` and then uses `XA`. Cannot run as stored. It builds the initial vector of the old
DAE version (3810 unknowns). *Status:* moved to `archive/root_DAE_and_scratch/`, not fixed.

## 3. `Mod/reversalEvent.m` uses an undefined global `nascar`
`nascar` is never defined anywhere, and `Mod/FEM_NL.m` does not call the function. The only place
that uses the event is the abandoned `Mod/DXI_algebric/FEM_NL.m`. *Status:* archived, not fixed.

## 4. `Loss/Furbo_Pre/` needs files that are not in its folder
`FEM1D.m` does `load('KFF.mat')` and `load('p.mat')` but neither file is in the folder (the sibling
`xiT0nullo/` has them). It also uses 19 Hz for the initial damping but 18 Hz inside its `myfun.m`.
*Status:* archived (`archive/loss_variants/Furbo_Pre/`), not fixed.

## 5. Force term of the linear transient case was lost (`Transient_linear/`)
`FEM1D.m` calls `fun_ode_1D(T,X,M,C,K,N)` (6 arguments) but the stored `fun_ode_1D.m` has 4 arguments and
no force term, so the load `FFF` set by `funforce.m` had no effect. The right-hand side
`M a = -K u - C v + F` with `F = FFF` is **inferred** (in `src/linear_beam_solve.m`), not recovered. The
stored parameters (37.5 N rectangular, 50 N incremental) also differ from the 200 N / 18 N s quoted in thesis
section 3.3.2, so it is unclear which exact runs the thesis figures came from.

## 6. `Brinson/Tension_compression_parametric.m` could not produce the stored `p.mat`
The file as stored had a single amplitude `eps_upp = 0.025` (the polynomial fit degenerates),
`eps_L = 0.0086 / -0.00746` (not the thesis Table 3.2 values 0.016 / -0.0146) and no `save`;
`Loss/p.mat` was saved by hand from some earlier run.
*Status:* the script now uses `eps_upp = 0.001:0.001:0.025` and `0.016 / -0.0146` (reconstructed from the
thesis) and writes `results/p_regen.mat`. The **original `data/p.mat` is shipped and used as
default**. The regenerated curve matches it within 0.001-0.003 for strains >= 0.3 % but is lower at
small strain (eta at 0.06 %: 0.119 versus 0.137 original; at 0.1 %: 0.131 versus 0.146), see
`docs/loss_factor_fit.png`. The original run therefore used a slightly different amplitude set (unknown).

## 7. Loss/KFF_Computation.m never saves; shipped `KFF.mat` files do not match it
It builds the initial stiffness matrix but has no `save`; the `KFF.mat` in the original folders was saved
by hand, and differs ~12 % from what the script computes. The main path does not need `KFF.mat`
(the solver builds the matrix), so the script was dropped (a copy is in `archive/loss_variants/xiT0nullo/`).

## 9. The "mid-span" displacement was the deflection of node 9, not of x = L/2 - FIXED

`FEM1D.m` read the mid-span displacement from `X(:,45)`, i.e. the deflection dof of node 9 (dof map of
`build_beam.m`: node k has w = dof 2k-3, theta = dof 2k-2). With N = 16 nodes (15 elements) the middle of the
beam is the middle of element 8 (between nodes 8 and 9), so node 9 sits at x = 533 mm and its peaks are
~0.56 % lower than at x = 500 mm for mode 1. `res.midspan` is now the deflection at exactly x = L/2, interpolated with
the cubic Hermite functions of the element from w and theta (theta = dw/dx, the same convention as the stiffness
and strain matrices) of nodes 8 and 9; if N-1 is even the middle node is used. The old value remains as
`res.midspan_node` / `res.x_midspan_node` / `res.col`. Base peaks, old -> new: 7.087 -> 7.128, 5.490 -> 5.521,
4.149 -> 4.172, 3.072 -> 3.090, 2.239 -> 2.252, 1.612 -> 1.621, 1.149 -> 1.155, 0.812 -> 0.817 mm. The state XX and
all physics are unchanged; only the read-out changed.

**Still as in the original (not changed):** the force of the `pulse` cases is applied to elements N/2 and N/2+1
(elements 8 and 9), i.e. symmetric around node 9 (x = 533 mm), not around L/2. The initial velocity field
(`initial_velocity.m`) is not affected by this. If a load exactly at mid-span is wanted, the force would have to be split between
nodes 8 and 9.

## 8. Smaller things
* `functxi.m` (archive): operator-precedence bug `xi_T(indp) .* 1-xi_sp(indp)`; never called.
* `Brinson/Prova_zone.m` calls `Zone` with the wrong number of arguments; `Zone_simpl.m` has a syntax
  error (`phase 'MD3n'`); `Prestrain_epsilon.m` has an index error. All three are in `archive/brinson_unfinished/`.
* `meshare.m` had the 30 x 20 mm section hard-coded and the fixed title `t = 0.8`; both are now
  arguments / computed (the base beam is 15 x 10 mm). It used the third-party `smooth2a.m` (Greg Reeves 2009, no licence);
  replaced by `src/smooth_mean2d.m` (same result to 4e-16).
* The base case has b = 15, t = 10 mm while most other cases use b = 30, t = 20 mm (see `cases/`).
* The reference frequency `f_ref` (10 / 14.5 / 18 / 19 / 21 Hz) had to be edited by hand in two files
  (`FEM1D.m` and `myfun.m`); it is now one parameter (`cfg.damping.f_ref`).
* The thesis abstract says the beam is integrated with `ode15s`; the final method (`Loss/`) uses `ode23t`
  (the abandoned DAE versions used `ode15s`). Only `prestrain_se` still uses `ode15s` in MATLAB.
* `ode15s` (SUNDIALS) in Octave fails at the first stress reversal of the pre-strain case ("error
  test failed"), therefore `prestrain_se` uses the `ode23t` stand-in under Octave.
* Original files were in Latin-1 / Windows-1252 and CRLF (14 of 114 `.m` files non-UTF-8: the degree sign in
  comments); converted to UTF-8 + LF.
