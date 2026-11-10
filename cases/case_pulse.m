function cfg = case_pulse()
% FORCED CASE  Rectangular force pulse on the SMA beam (thesis 3.3.2, Fig. 3.23).
%   Beam 30 x 20 mm. 200 N in total (2 x 100 N on the two mid-span elements), ramp 0 -> F in
%   1e-4 s, held until t = 0.02 s, then released; starts at rest. f_ref = 18 Hz, loss factor at 0.001.
%   Original: Loss/Transient/FEM1D.m (+ its own myfun.m / fun_ode_1D.m).
%   SLOW: 0.16 s of the 0.8 s took ~33 min in Octave (see README: runtimes).
cfg = default_config();
cfg.name = 'pulse';
cfg.description = 'Rectangular force pulse 200 N / 0.02 s on the 30x20 mm SMA beam (Loss/Transient)';
cfg.geom.b = 30;  cfg.geom.thick = 20;
cfg.damping.f_ref = 18;
cfg.damping.x_ref = 0.001;
cfg.init.type = 'zero';
cfg.load.type = 'pulse';
cfg.load.F_element = -100;     % [N] per element -> 200 N total
cfg.load.t_rise = 1e-4;
cfg.load.t_end  = 0.02;
cfg.solver.InitialStep = 5e-6;
cfg.post.shift_to_first_max = false;
end
