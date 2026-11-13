function cfg = case_composite()
% COMPOSITE  GFRP / SMA / GFRP sandwich beam (thesis chapter 4.1, configuration 1).
%   SMA core 15 x 10 mm (nonlinear) + two GFRP skins (linear, total thickness 12 mm).
%   Reference frequency 14.5 Hz, loss factor taken at strain 0.0007.
cfg = default_config();
cfg.name = 'composite';
cfg.description = 'GFRP/SMA/GFRP sandwich, free vibration';
cfg.damping.f_ref = 14.5;
cfg.damping.x_ref = 0.0007;
cfg.init.scale = 0.663;
cfg.hle.N = cfg.geom.N;  cfg.hle.L = cfg.geom.L;  cfg.hle.b = 15;
cfg.hle.rho = 1.88e-9;  cfg.hle.E = 46e3;       % GFRP
cfg.hle.thick = 12;     cfg.hle.thick2 = 10;     % total / SMA-core thickness [mm]
cfg.hle.rho2 = 6.1e-09; cfg.hle.E2 = 26e3;       % SMA (linear, only for the Rayleigh coefficients)
cfg.hle.zeta = 0.013;                            % Rayleigh ratio on modes 1 and 2
end
