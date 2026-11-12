function cfg = default_config()
% DEFAULT_CONFIG  every parameter a simulation can take, with the thesis base-case values.
%
% A "case" is just a function cases/case_<name>.m that starts from this struct and changes
% what is different. run_case('<name>') reads the struct and calls the solver, so a future
% interactive menu only has to fill in the fields marked  [MENU]  and call
% run_case(cfg)  with the struct.
%
% Units: mm, s, MPa, ton/mm^3 (so that N, mm, MPa, s are consistent), Hz.

cfg.name        = 'default';
cfg.description = '';
cfg.model       = 'two_phase';   % 'two_phase' (ferro-elastic, xi_s+, xi_s-, xi_T) | 'one_phase_SE' (super-elastic) | 'linear'

% ---- geometry ---------------------------------------------------------------------------
cfg.geom.N     = 16;      % [MENU] nodes (N-1 = 15 elements). Only N=16 is verified; other N untested
cfg.geom.L     = 1000;    % beam length [mm]
cfg.geom.b     = 15;      % width [mm]
cfg.geom.thick = 10;      % thickness [mm]
cfg.geom.n_csi = 5;       % integration points along one element
cfg.geom.n_y   = 10;      % integration points through the thickness  (N-1)*5*10 = 750 points

% ---- SMA material (Cu-Al-Zn-Mn, thesis Table 3.2) ---------------------------------------------------
cfg.mat.Ea     = 34.9e3;  % austenite Young's modulus [MPa]
cfg.mat.Emp    = 26e3;    % detwinned martensite, tension
cfg.mat.Emm    = 33e3;    % detwinned martensite, compression
cfg.mat.EmT    = 26e3;    % twinned martensite
cfg.mat.ssp    = 0;       % critical stresses [MPa]: start / finish of tensile detwinning
cfg.mat.sfp    = 188;
cfg.mat.ssm    = 0;       % start / finish of compressive detwinning
cfg.mat.sfm    = -220;
cfg.mat.eps_Lp = 0.016;   % max residual strain, tension
cfg.mat.eps_Lm = -0.0146; % max residual strain, compression
cfg.mat.rho    = 6.1e-09; % density [ton/mm^3]
cfg.mat.E_K    = 26e3;    % modulus used for the first (linear) stiffness matrix (= EmT)

% ---- loss factor curve ----------------------------------------------------------------------------
cfg.p_file     = 'p.mat'; % in data/: polynomial eta(strain) from brinson/Tension_compression_parametric.m
cfg.damping.f_ref = 10;   % reference frequency in C = eta/(2*pi*f_ref)*K  [Hz]   [MENU: frequency]
cfg.damping.x_ref = 0.0006; % strain at which the initial damping is evaluated

% ---- initial state ------------------------------------------------------------------------------------
cfg.init.type  = 'mode1_velocity';  % half-sine initial velocity field | 'zero'
cfg.init.scale = 0.46;               % amplitude factor of the initial velocity field

% ---- external load  (type 'none' = free vibration) -----------------------------------------------------
cfg.load.type      = 'none';   % [MENU] 'none' | 'pulse' (ramp up, hold, release)
cfg.load.F_element = 0;        % [MENU] nodal force [N] on each of the two mid-span elements (total = 2x)
cfg.load.t_rise    = 1e-4;     % ramp-up time [s]
cfg.load.t_end     = 0.02;     % force released at this time [s]
cfg.load.dt_restart = 1e-5;    % solver restart window after each change of slope [s]

% ---- time ---------------------------------------------------------------------------------------------------
cfg.time.Tf    = 0.8;      % [MENU] duration [s]

% ---- ODE solver -------------------------------------------------------------------------------------------
cfg.solver.name         = 'ode23t';   % MATLAB
cfg.solver.name_octave  = 'ode23t';   % Octave: compat/octave/ode23t.m forwards to ode23s
cfg.solver.InitialStep  = [];         % only for the first call
cfg.solver.RelTol       = [];
cfg.solver.AbsTol       = [];
cfg.solver.max_restarts = 1e6;

% ---- numerics -------------------------------------------------------------------------------------------------
cfg.numerics.newton_maxit = 100;    % cap on the Newton iterations per step (the loop runs while ANY point has not converged, || - see KNOWN_ISSUES.md #1)
cfg.numerics.verbose      = false;  % true = original per-step printout

% ---- post-processing --------------------------------------------------------------------------------------------
cfg.post.shift_to_first_max = true; % time axis starts at the first mid-span maximum (as in the original)
cfg.post.keep_fields        = false; % keep the strain / martensite fields of the LAST step (used by plot_fields)
cfg.post.keep_field_history = false; % with keep_fields: keep them for EVERY step (large and slow, O(n^2))

% ---- performance switches (see PERF_NOTES.md). 1 = fast path, 0 = original code path (identical physics) ----
cfg.perf.vec      = 1;   % vectorised element integration/assembly (replaces the trapz loops) and strain B*Extr
cfg.perf.sets     = 1;   % logical masks instead of setdiff/union/ismember in the OutputFcn
cfg.perf.hist     = 0;   % 1 = record tt, yy and the per-step strain/martensite fields (never read back; slow O(n^2))
cfg.perf.fast_rhs = 1;   % constant inv(M) instead of M\ at every right-hand-side call
cfg.perf.jac      = 1;   % analytic Jacobian [-Minv*C -Minv*K; I 0] handed to the solver

% ---- extra elastic layers (composite only) --------------------------------------------------------
cfg.hle = [];
end
