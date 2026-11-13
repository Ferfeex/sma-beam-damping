function cfg = case_base()
% BASE CASE  Free vibration of the SMA beam (thesis chapter 3.3.2; the ring-down is the one of thesis Fig. 4.3).
%   Cu-Al-Zn-Mn beam 1000 x 15 x 10 mm, simply supported, ferro-elastic (T below Mf),
%   released with an initial half-sine velocity field. It rings down because every stress
%   reversal dissipates energy (loss factor eta ~ 0.14, C = eta/(2*pi*10 Hz) K).
%   Original: Loss/FEM1D.m (+ myfun.m, fun_ode_1D.m).  Expected peaks of the mid-span
%   displacement at x = L/2 (Hermite) after the time shift: 7.13 (t=0), 5.52, 4.17, 3.09, 2.25, 1.62 mm
%   (the node-9 value (x = 533 mm) used by earlier versions is 0.56 % lower: 7.09, 5.49, 4.15, 3.07, 2.24, 1.61).
cfg = default_config();
cfg.name = 'base';
cfg.description = 'Free vibration, SMA beam 15x10 mm, ferro-elastic, loss-factor damping (Loss/FEM1D.m)';
end
