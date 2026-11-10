function cfg = case_base()
% BASE CASE  Free vibration of the SMA beam (thesis chapter 3.3.2, Fig. 3.x / 4.3).
%   Cu-Al-Zn-Mn beam 1000 x 15 x 10 mm, simply supported, ferro-elastic (T below Mf),
%   released with an initial half-sine velocity field. It rings down because every stress
%   reversal dissipates energy (loss factor eta ~ 0.14, C = eta/(2*pi*10 Hz) K).
%   Original: Loss/FEM1D.m (+ myfun.m, fun_ode_1D.m).  Expected peaks of the mid-span
%   displacement after the time shift: 7.09 (t=0), 5.49, 4.16, 3.08, 2.24, 1.61 mm.
cfg = default_config();
cfg.name = 'base';
cfg.description = 'Free vibration, SMA beam 15x10 mm, ferro-elastic, loss-factor damping (Loss/FEM1D.m)';
end
