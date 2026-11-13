function cfg = case_hle_linear()
% LINEAR REFERENCE ("HLE")  Same beam as case_base but with a linear elastic material of
%   E = 26 GPa (= the martensite modulus) and Rayleigh damping (1.3 % on modes 1 and 2).
%   Used to show what the SMA damping adds (E = 26e3, rho = 6.1e-9).
cfg = default_config();
cfg.name = 'hle_linear';
cfg.description = 'Linear reference beam, E=26 GPa, Rayleigh 1.3 %, same initial velocity as the base case';
cfg.model = 'linear';
cfg.mat.E = 26e3;
cfg.damping.zeta = 0.013;
cfg.solver.name = 'ode23t';
cfg.post.shift_to_first_max = true;
end
