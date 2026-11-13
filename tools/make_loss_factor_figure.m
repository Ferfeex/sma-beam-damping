% MAKE_LOSS_FACTOR_FIGURE  docs/loss_factor_fit.png: loss factor eta versus reversal strain,
% computed by brinson/Tension_compression_parametric.m (dots), its 4th-degree polynomial fit
% (regenerated p) and the original data/p.mat that the simulations use by default.
root = fileparts(fileparts(mfilename('fullpath')));
set(0,'defaultfigurevisible','off');
setenv('SMA_ROOT',root);        % the script below does 'clear all'
cd(fullfile(root,'brinson'));
run('Tension_compression_parametric.m');          % creates p, eu, Loss_factor; writes results/p_regen.mat
root = getenv('SMA_ROOT');
orig = load(fullfile(root,'data','p.mat'));
f = figure('visible','off');
xx = linspace(0,0.025,300);
plot(eu*100, Loss_factor(2:end), 'ko'); hold on
plot(xx*100, polyval(p,xx), 'LineWidth', 2);
plot(xx*100, polyval(orig.p,xx), '--', 'LineWidth', 2);
grid on; xlabel('Upper reversal strain [%]'); ylabel('Loss factor \eta');
title('Loss factor of Cu-Al-Zn-Mn from hysteresis cycles (Brinson model)');
legend('cycle-by-cycle \eta = E_{diss}/(2\pi U_{max})','4th-degree fit, regenerated now','original data/p.mat (used by default)','Location','southeast');
print(f,'-dpng','-r110',fullfile(root,'docs','loss_factor_fit.png'));
