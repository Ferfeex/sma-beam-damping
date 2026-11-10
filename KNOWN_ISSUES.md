# Known issues of the original (2020) code

Found while reading, porting and cleaning the code. "Status" says what was done about each one.
Nothing here changes the default results of the base case.

## 1. Suspected `&&` slip in the Newton loop (`src/myfun.m`, `src/myfun_SE.m`)

```matlab
while ~isempty(indp) && ~isempty(indm)      % original
```

`indp` = integration points that are in the tension transformation zone, `indm` = points in the
compression zone. With `&&` the Newton iteration **runs only while both sets are non-empty**, and it
**stops as soon as either set has converged**, even if the other one has not. If only one set is
active at a time step (for instance only tension zones), the martensite fractions of those points are
not updated at all in that step. An `||` ("while any point still has to converge") looks like the
intended condition.

**Status: reported, NOT changed.** The default is the original `&&` (so all results are the original
ones); the `||` version is an opt-in flag, `run_case('base','newton_or',true)` (with an iteration cap of
100 so that a non-converging point cannot hang the loop).

**Before / after, full base case (0.8 s, same Octave, same integrator):**

| | original `&&` | `\|\|` fix |
|---|---|---|
| solver steps | 4320 | 4883 |
| run time | 370 s | 431 s |
| mid-span peaks [mm] | 7.087, 5.494, 4.155, 3.076, 2.239, 1.611, 1.147, 0.811 | 7.087, 5.490, 4.149, 3.075, 2.241, 1.612, 1.150, 0.814 |
| peak times [s] | 0, .1105, .2203, .3287, .4371, .5441, .6510, .7562 | 0, .1104, .2200, .3288, .4367, .5428, .6499, .7570 |
| max difference of the mid-span curve | | 0.058 mm (0.8 % of the 7.09 mm amplitude), rms 0.013 mm |

(Run times were measured under different machine load, so they are indicative only.)

So the effect on the result is tiny (peak heights differ by <= 0.006 mm), of the same size as the
differences that a different step sequence causes anyway, but `||` needs ~13 % more steps. The
conclusion for the thesis does not change; it is Federico's call whether to adopt `||`.
(The composite and the other cases were not re-run with `||`.)

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
thesis) and writes `results/p_regen.mat`. The author's **original `data/p.mat` is shipped and used as
default**. The regenerated curve matches it within 0.001-0.003 for strains >= 0.3 % but is lower at
small strain (eta at 0.06 %: 0.119 versus 0.137 original; at 0.1 %: 0.131 versus 0.146), see
`docs/loss_factor_fit.png`. The original run therefore used a slightly different amplitude set (unknown).

## 7. Loss/KFF_Computation.m never saves; shipped `KFF.mat` files do not match it
It builds the initial stiffness matrix but has no `save`; the `KFF.mat` in the original folders was saved
by hand, and differs ~12 % from what the script computes. The main path does not need `KFF.mat`
(the solver builds the matrix), so the script was dropped (a copy is in `archive/loss_variants/xiT0nullo/`).

## 8. Smaller things
* `functxi.m` (archive): operator-precedence bug `xi_T(indp) .* 1-xi_sp(indp)`; never called.
* `Brinson/Prova_zone.m` calls `Zone` with the wrong number of arguments; `Zone_simpl.m` has a syntax
  error (`phase 'MD3n'`); `Prestrain_epsilon.m` has an index error. All three are in `archive/brinson_unfinished/`.
* `meshare.m` had the 30 x 20 mm section hard-coded and the fixed title `t = 0.8`; both are now
  arguments / computed (the base beam is 15 x 10 mm).
* The base case has b = 15, t = 10 mm while most other cases use b = 30, t = 20 mm (see `cases/`).
* The reference frequency `f_ref` (10 / 14.5 / 18 / 19 / 21 Hz) had to be edited by hand in two files
  (`FEM1D.m` and `myfun.m`); it is now one parameter (`cfg.damping.f_ref`).
* The thesis abstract says the beam is integrated with `ode15s`; the final method (`Loss/`) uses `ode23t`
  (the abandoned DAE versions used `ode15s`). Only `prestrain_se` still uses `ode15s` in MATLAB.
* `ode15s` (SUNDIALS) in Octave fails at the first stress reversal of the pre-strain case ("error
  test failed"), therefore `prestrain_se` uses the `ode23t` stand-in under Octave.
* Original files were in Latin-1 / Windows-1252 and CRLF (14 of 114 `.m` files non-UTF-8: the degree sign in
  comments); converted to UTF-8 + LF.
