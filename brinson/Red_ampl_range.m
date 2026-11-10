% Prove reduction amplitude range Brinson old

clear
clc
close all

% Material properties
E_a = 34.9e3; %[MPa]
E_m = 26e3; %[MPa]
Teta = 0.55; %[MPa/°C]
Mf = -55; %[°C]
Ms = -35; %[°C]
As = -30; %[°C]
Af = -15; %[°C]
C_M = 1; %[MPa/°C]
C_A = 3; %[MPa/°C]
sigma_s_cr = 70; %[MPa]
sigma_f_cr = 188; %[MPa]
eps_L = 0.016;
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);

T = 0;

sigma = linspace(0, C_M * (T - Ms) + sigma_s_cr, 1000)';
eps = sigma / E_a;
csi = zeros(1000,1);

sigma2 = linspace(C_M * (T - Ms) + sigma_s_cr, C_M * (T - Ms) + sigma_f_cr, 1000)';
csi2 = 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma2 - sigma_f_cr - C_M * (T - Ms))) + 1/2;
E2 = E_a + csi2 * (E_m - E_a);
eps2 = eps_L * csi2 + sigma2 ./ E2;

sigma3 = linspace(C_M * (T - Ms) + sigma_f_cr, C_A * (T - As), 1000)';
eps3 = eps_L + sigma3 / E2(end);
csi3 = ones(1000,1);

sigma4 = linspace(C_A * (T - As), C_A * (T - Af), 1000)';
csi4 = 1/2 * (cos(a_A * (T - As - sigma4 / C_A)) + 1);
E4 = E_a + csi4 * (E_m - E_a);
eps4 = eps_L * csi4 + sigma4 ./E4;


sigma5 = linspace(C_A * (T - Af), C_M * (T - Ms) + sigma_s_cr, 1000)';
eps5 = sigma5 / E_a;
csi5 = zeros(1000,1);

sigma6 = linspace(C_M * (T - Ms) + sigma_s_cr, C_M * (T - Ms) + sigma_f_cr - 30, 1000)';
csi6 = 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma6 - sigma_f_cr - C_M * (T - Ms))) + 1/2;
E6 = E_a + csi6 * (E_m - E_a);
eps6 = eps_L * csi6 + sigma6 ./ E6;

sigma7 = linspace(sigma6(end), C_A * (T - As), 1000)';
eps7= eps_L * csi6(end) + sigma7 / E6(end);
csi7 = ones(1000,1) * csi6(end);

sigma8 = linspace(C_A * (T - As), C_A * (T - Af) + 80, 1000)';
csi8 = csi6(end)/2 * (cos(a_A * (T - As - sigma8 / C_A)) + 1);
E8 = E_a + csi8 * (E_m - E_a);
eps8 = eps_L * csi8 + sigma8 ./E8;

sigma9 = linspace(sigma8(end),C_M * (T - Ms) + sigma_s_cr, 1000)';
eps9 = eps_L * csi8(end) + sigma9 ./E8(end);
csi9 = ones(1000,1) * csi8(end);

sigma10 = linspace(C_M * (T - Ms) + sigma_s_cr, C_M * (T - Ms) + sigma_f_cr, 1000)';
csi10 = (1 - csi9(end))/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma10 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi9(end))/2;
E10 = E_a + csi10 * (E_m - E_a);
eps10 = eps_L * csi10 + sigma10 ./ E10;

% sigma11 = linspace(sigma10(end), C_A * (T - As), 1000)';
% eps11= eps_L * csi10(end) + sigma11 / E10(end);
% csi11 = ones(1000,1) * csi10(end);

% sigma12 = linspace(sigma11(end), C_A * (T - As), 1000)';
% csi12 = csi10(end)/2 * (cos(a_A * (T - As - sigma12 / C_A)) + 1);
% E12 = E_a + csi12 * (E_m - E_a);
% eps12 = eps_L * csi12 + sigma12 ./E12;

sigma = [sigma; sigma2; sigma3; sigma4; sigma5; sigma6; sigma7; sigma8; sigma9; sigma10];
eps = [eps; eps2; eps3; eps4; eps5; eps6; eps7; eps8; eps9; eps10];
csi = [csi; csi2; csi3; csi4; csi5; csi6; csi7; csi8; csi9; csi10];

figure
plot(sigma, csi)
grid on
xlabel('\sigma [MPa]')
ylabel('\xi')

figure
plot(eps, sigma)
grid on
xlabel('\epsilon')
ylabel('\sigma [MPa]')
