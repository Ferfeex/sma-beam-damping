function cfg = case_prestrain_se()
% PRE-STRAIN / SUPER-ELASTIC CASE  (thesis 3.3.3, Figs 3.25-3.30)
%   SMA beam 30 x 20 mm above Af (T0 = -10 degC in the thesis phase-diagram convention),
%   pre-stressed with sigma0 = 60 MPa, one martensite fraction (xi_s+), f_ref = 21 Hz.
%   The loss-factor curve here is a function of STRESS (data/p_SE.mat, from brinson/PrestrainSE.m).
%   MATLAB uses ode15s; Octave's ode15s (SUNDIALS) fails at the first stress reversal, so in
%   Octave the ode23t substitute is used (see cfg.solver.name_octave).
cfg = default_config();
cfg.name = 'prestrain_se';
cfg.description = 'Super-elastic beam with 60 MPa pre-stress';
cfg.model = 'one_phase_SE';
cfg.geom.b = 30;  cfg.geom.thick = 20;
cfg.mat.ssp = 35;
cfg.mat.Mf = -55; cfg.mat.Ms = -35; cfg.mat.As = -30; cfg.mat.Af = -15;
cfg.mat.CMp = 1;  cfg.mat.CAp = 3;
cfg.p_file = 'p_SE.mat';
cfg.damping.f_ref = 21;
cfg.damping.x_ref = 80;          % polyval(p, 80): initial loss factor evaluated at 80 MPa
cfg.damping.LF_offset = 0.015;
cfg.damping.sigma_rev = 120;
cfg.init.scale = 1;
cfg.init.Temp0 = -10;
cfg.init.sigma0 = 60;
cfg.solver.name = 'ode15s';
cfg.solver.name_octave = 'ode23t';
cfg.solver.RelTol = 1e-2;
cfg.solver.AbsTol = 1e-5;
end
