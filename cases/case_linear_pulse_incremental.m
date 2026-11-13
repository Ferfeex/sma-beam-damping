function cfg = case_linear_pulse_incremental()
% LINEAR GFRP BEAM, INCREMENTAL (ramp) PULSE  (thesis 3.2 / Fig. 3.21)
%   GFRP-like linear beam (E=46 GPa, rho=1.88e-9), Rayleigh 1.3 %. The force ramps linearly
%   from 0 to 2*25 N at t = 0.033 s and is then released. The force history is
%   reconstructed from the thesis figure (see src/linear_beam_solve.m).
cfg = default_config();
cfg.name = 'linear_pulse_incremental';
cfg.description = 'Linear GFRP beam, ramp load 50 N to 0.033 s then released';
cfg.model = 'linear';
cfg.mat.E = 46e3;  cfg.mat.rho = 1.88e-09;
cfg.damping.zeta = 0.013;
cfg.init.type = 'zero';
cfg.load.type = 'ramp';
cfg.load.F_element = -25;
cfg.load.t_rise = 0.033; cfg.load.t_end = 0.033; cfg.load.dt_restart = 1e-4;
cfg.solver.InitialStep = 1e-4;
cfg.post.shift_to_first_max = false;
end
