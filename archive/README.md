# archive/ - older attempts and side paths (NOT maintained)

Everything in here is kept for history and because it shows how the project evolved, but it
is **not part of the working code** and is not expected to run. Only the `.m` sources were
kept: the large data/result files (`.mat` of 5-170 MB, ~1 GB in total) and the `.fig` plots
stayed on the original Drive folder. Most scripts below therefore also depend on data files
that are not here.

Folder names are the original ones, so you can map them to `docs/CODE_MAP_original_analysis.md`.

| Folder | What it was | Status |
|---|---|---|
| `ODE_NEWT_first_nonlinear_approach/` | **First nonlinear (SMA) beam** (May-Jun 2020): a Newton / `fsolve` solve of the martensite fractions *inside* the ODE right-hand side. Correct idea but very slow; superseded by the OutputFcn method in `src/`. `Variabili_1000N/` post-processes a stored run (`track.m`, `evolution.m`). Some files use a different parameter set (Emp 130/140 GPa ...). | needs `matlab.mat` etc. (not included) |
| `Mod_DAE_ode15s/` | **Intermediate version** (Jun 2020): the nonlinear problem written as a DAE for `ode15s` (3810 unknowns: dynamics + one algebraic constraint per integration point). `DXI_algebric/` is a variant with 6060 unknowns (martensite-fraction *rates* as unknowns, `odeset('Events',@reversalEvent)`); it appears abandoned. | needs `Sparsity.mat`, `IN_COND.mat` ... (not included); `reversalEvent.m` uses an undefined global `nascar` |
| `root_DAE_and_scratch/` | Latest copy of the DAE version (`FEM_NL.m`, `fun_ode_NL.m`, `fun_kff.m`, `fun_xi.m`, `fun_eps.m`, `inizcon.m`) plus loose scratch files (`functxi.m`, `myfun2.m`, `funforce.m`, `reversalEvent.m`, `meshare.m` ...). | `Condizioni_iniziali_unforced.m` is broken (uses undefined `XA`); `functxi.m` has an operator-precedence bug (`xi_T .* 1-xi_sp`) and is never called |
| `FEM2D/` | Early experiment: 2D beam element (axial + bending) with a harmonic force. Not connected to the rest and not used in the thesis text. | standalone, untested |
| `loss_variants/Furbo_Pre/` | Same method as the base case with b=30, t=20, f=19 Hz, initial fractions xi_s+0 = xi_s-0 = 0.5, xi_T0 = 0 ("furbo" = Italian "cunning"). Uses a stress-offset formulation (`sigma0`, `epsilon0`) in `myfun.m`. | needs `KFF.mat` and `p.mat` that were not in the folder |
| `loss_variants/xiT0nullo/` | Same idea with xi_T0 null; only differs from `Furbo_Pre` in Emm and E. Probably the "xi_T0 variation" of thesis Fig. 3.20 (unclear). Has its own `KFF.mat` and `p.mat` (small, kept). | probably runnable from inside the folder with the Octave `ode23t` stand-in (not tested) |
| `loss_variants/meshare_variants/` | Field-plot variants for the SE and pulse cases (the base `meshare.m` is in `src/`). | |
| `brinson_unfinished/` | Material-level helpers that do not run: `Prova_zone.m` (calls `Zone` with the wrong arguments), `Zone.m`/`Zone_simpl.m` (`Zone_simpl` has a syntax error), `Prestrain_epsilon.m` (index error), `Newton.m`/`Newton_sempl.m` (only referenced from comments). | broken |
| `original_case_scripts/` | The **unmodified originals** (only converted to UTF-8 + LF) of the scripts that were turned into `cases/*.m`: `Loss/FEM1D.m` (base), Composite, SE_pre-strain, Transient, Transient_linear (incremental + rectangular), root `FEM1D.m` (linear), and their `myfun`/`fun_ode_1D` variants. | used by `tests/run_tests.m` to prove the refactor changed nothing; diff them with `src/` to see exactly what changed |

## Duplicates removed

The original Drive folder held 114 `.m` files of which only 94 were distinct (20 byte-identical
copies: e.g. `smooth2a.m` x7, `Loss/*/fun_ode_1D.m` x5, `meshare.m` x4 ...). The copies were dropped,
one copy of each was kept (in `src/`, `third_party/` or here). `CLEANUP_REPORT.md` lists the mapping.

Note: some archived scripts (`meshare.m`, `track.m` ...) call `smooth2a`, a third-party function that is no longer in
the repository (replaced in the main code by `src/smooth_mean2d.m`); they were not runnable anyway.
