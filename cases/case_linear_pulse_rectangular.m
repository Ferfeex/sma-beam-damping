function cfg = case_linear_pulse_rectangular()
% LINEAR GFRP BEAM, RECTANGULAR PULSE  (thesis 3.2 / Fig. 3.23)
%   As case_linear_pulse_incremental but a rectangular pulse of 2*18.75 N, ramp 1e-4 s,
%   released at t = 0.022 s.
cfg = default_config();
cfg.name = 'linear_pulse_rectangular';
cfg.description = 'Linear GFRP beam, rectangular pulse 37.5 N until 0.022 s';
cfg.model = 'linear';
cfg.mat.E = 46e3;  cfg.mat.rho = 1.88e-09;
cfg.damping.zeta = 0.013;
cfg.init.type = 'zero';
cfg.load.type = 'pulse';
cfg.load.F_element = -18.75;
cfg.load.t_rise = 1e-4; cfg.load.t_end = 0.022; cfg.load.dt_restart = 1e-5;
cfg.solver.InitialStep = 1e-5;
cfg.post.shift_to_first_max = false;
end
