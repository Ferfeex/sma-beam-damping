clear
clc
close all
% Cu-Al-Zn-Mn
E_a = 34.9e3; %[MPa]
E_m_p = 26e3; %[MPa]
E_m_n = 33e3;
E_m_T = 26e3;
% E_m_T+ = E_m_T- --> Ricorda che twinned negativa = twinned positiva
Teta = 0.55; %[MPa/°C]
Mf = -55; %[°C]
Ms = -35; %[°C]
As = -30; %[°C]
Af = -15; %[°C]
C_M_p = 1; %[MPa/°C]
C_A_p = 3; %[MPa/°C]
C_M_n = 0.2;
C_A_n = 3;
sigma_s_cr_p = 70; %[MPa]
sigma_f_cr_p = 188; %[MPa]
sigma_s_cr_n = -0;
sigma_f_cr_n = -220;
eps_L_p = 0.016;
eps_L_n = -0.0146;
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);
T = -60;
%% Tensione-compressione
E = E_m_T;
csi_S_p = 0;
csi_S_n = 0;
csi_T = 1;
% Linear Tension
sigma = linspace(0, sigma_s_cr_p, 1000)';
eps = sigma / E(end) + eps_L_p * csi_S_p(end) + eps_L_n * csi_S_n(end);
% Reorientation T
sigma2 = linspace(sigma(end), 120, 1000)';
Y2 = cos(pi * (sigma2 - sigma_s_cr_p) / (sigma_f_cr_p - sigma_s_cr_p));
csi_S_p2 = (1 - Y2) / 2 + csi_S_p(end) * (1 + Y2) / 2;
csi_T2 = csi_T(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
csi_S_n2 = csi_S_n(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
E2 = E_a + csi_S_p2 * (E_m_p - E_a) + csi_S_n2 * (E_m_n - E_a) + csi_T2 * (E_m_T - E_a);
eps2 = sigma2 ./ E2 + eps_L_p * csi_S_p2 + eps_L_n * csi_S_n2;
eps_down = -eps2(end);
% Linear Unloading
sigma3 = linspace(sigma2(end), sigma_s_cr_p, 1000)';
eps3 = sigma3 / E2(end) + eps_L_p * csi_S_p2(end) + eps_L_n * csi_S_n2(end);
% Self-Accommodation
sigma4 = linspace(sigma3(end), sigma_s_cr_n, 1000)';
Y4 = cos(pi * (sigma4 - sigma_s_cr_p) / (sigma_f_cr_n - sigma_s_cr_p));
csi_T4 = (1 - Y4) / 2 + csi_T2(end) * (1 + Y4) / 2;
csi_S_n4 = csi_S_n2(end);
csi_S_p4 = csi_S_p2(end) * (1 - csi_T4 - csi_S_n4) / (1 - csi_T2(end) - csi_S_n4);
E4 = E_a + csi_S_p4 * (E_m_p - E_a) + csi_S_n4 * (E_m_n - E_a) + csi_T4 * (E_m_T - E_a);
eps4 = sigma4 ./ E4 + eps_L_p * csi_S_p4 + eps_L_n * csi_S_n4;
% Reorientation C
sigma5 = linspace(sigma4(end), sigma_f_cr_n, 1000)';
Y5 = cos(pi * (sigma5 - sigma_s_cr_n) / (sigma_f_cr_n - sigma_s_cr_n));
csi_S_n5 = (1 - Y5) / 2 + csi_S_n4(end) * (1 + Y5) / 2;
csi_S_p5 = csi_S_p4(end) * (1 - csi_S_n5) / (1 - csi_S_n4(end));
csi_T5 = csi_T4(end) * (1 - csi_S_n5) / (1 - csi_S_n4(end));
E5 = E_a + csi_S_p5 * (E_m_p - E_a) + csi_S_n5 * (E_m_n - E_a) + csi_T5 * (E_m_T - E_a);
eps5 = sigma5 ./ E5 + eps_L_p * csi_S_p5 + eps_L_n * csi_S_n5;
eps5 = eps5(eps5 >= eps_down);
sigma5 = sigma5(1:length(eps5));
csi_S_p5 = csi_S_p5(1:length(eps5));
csi_S_n5 = csi_S_n5(1:length(eps5));
csi_T5 = csi_T5(1:length(eps5));
E5 = E5(1:length(eps5));
% Linear C Unloading
sigma6 = linspace(sigma5(end), sigma_s_cr_n, 1000)';
eps6 = sigma6 / E5(end) + eps_L_p * csi_S_p5(end) + eps_L_n * csi_S_n5(end);
%Self-Accommodation
sigma7 = linspace(sigma6(end), sigma_s_cr_p, 1000)';
Y7 = cos(pi * (sigma7 - sigma_s_cr_n) / (sigma_f_cr_p - sigma_s_cr_n));
csi_T7 = (1 - Y7) / 2 + csi_T5(end) * (1 + Y7) / 2;
csi_S_p7 = csi_S_p5(end);
csi_S_n7 = csi_S_n5(end) * (1 - csi_T7 - csi_S_p7) / (1 - csi_T5(end) - csi_S_p7);
E7 = E_a + csi_S_p7 * (E_m_p - E_a) + csi_S_n7 * (E_m_n - E_a) + csi_T7 * (E_m_T - E_a);
eps7 = sigma7 ./ E7 + eps_L_p * csi_S_p7 + eps_L_n * csi_S_n7;
% Reorientation T 2
sigma8 = linspace(sigma7(end), 120, 1000)';
Y8 = cos(pi * (sigma8 - sigma_s_cr_p) / (sigma_f_cr_p - sigma_s_cr_p));
csi_S_p8 = (1 - Y8) / 2 + csi_S_p7(end) * (1 + Y8) / 2;
csi_T8 = csi_T7(end) * (1 - csi_S_p8) / (1 - csi_S_p7(end));
csi_S_n8 = csi_S_n7(end) * (1 - csi_S_p8) / (1 - csi_S_p7(end));
E8 = E_a + csi_S_p8 * (E_m_p - E_a) + csi_S_n8 * (E_m_n - E_a) + csi_T8 * (E_m_T - E_a);
eps8 = sigma8 ./ E8 + eps_L_p * csi_S_p8 + eps_L_n * csi_S_n8;

sigma_plot = [sigma;sigma2;sigma3;sigma4;sigma5;sigma6;sigma7;sigma8];
eps_plot = [eps;eps2;eps3;eps4;eps5;eps6;eps7;eps8];

%% Loading-Unloading
E = E_m_T;
csi_S_p = 0;
csi_S_n = 0;
csi_T = 1;
% Linear Tension
sigma = linspace(0, sigma_s_cr_p, 1000)';
eps = sigma / E(end) + eps_L_p * csi_S_p(end) + eps_L_n * csi_S_n(end);
% Reorientation T
sigma2 = linspace(sigma(end), 120, 1000)';
Y2 = cos(pi * (sigma2 - sigma_s_cr_p) / (sigma_f_cr_p - sigma_s_cr_p));
csi_S_p2 = (1 - Y2) / 2 + csi_S_p(end) * (1 + Y2) / 2;
csi_T2 = csi_T(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
csi_S_n2 = csi_S_n(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
E2 = E_a + csi_S_p2 * (E_m_p - E_a) + csi_S_n2 * (E_m_n - E_a) + csi_T2 * (E_m_T - E_a);
eps2 = sigma2 ./ E2 + eps_L_p * csi_S_p2 + eps_L_n * csi_S_n2;
eps_down = -eps2(end);
% Linear Unloading
sigma3 = linspace(sigma2(end), sigma_s_cr_p, 1000)';
eps3 = sigma3 / E2(end) + eps_L_p * csi_S_p2(end) + eps_L_n * csi_S_n2(end);
% Self-Accommodation
sigma4 = linspace(sigma3(end), 0, 1000)';
Y4 = cos(pi * (sigma4 - sigma_s_cr_p) / (sigma_f_cr_n - sigma_s_cr_p));
csi_T4 = (1 - Y4) / 2 + csi_T2(end) * (1 + Y4) / 2;
csi_S_n4 = csi_S_n2(end);
csi_S_p4 = csi_S_p2(end) * (1 - csi_T4 - csi_S_n4) / (1 - csi_T2(end) - csi_S_n4);
E4 = E_a + csi_S_p4 * (E_m_p - E_a) + csi_S_n4 * (E_m_n - E_a) + csi_T4 * (E_m_T - E_a);
eps4 = sigma4 ./ E4 + eps_L_p * csi_S_p4 + eps_L_n * csi_S_n4;
% Reverse Self-Accommodation
sigma5 = linspace(sigma4(end), sigma_s_cr_p+(sigma2(end)-sigma_s_cr_p)/1.5, 1000)';
eps5 = sigma5 / E4(end) + eps_L_p * csi_S_p4(end) + eps_L_n * csi_S_n4(end);

sigma6 = linspace(sigma5(end), sigma_f_cr_p, 1000)';
Y6 = cos(pi * (sigma6 - sigma6(1)) / (sigma_f_cr_p - sigma6(1)));
csi_S_p6 = (1 - Y6) / 2 + csi_S_p4(end) * (1 + Y6) / 2;
csi_T6 = csi_T4(end) * (1 - csi_S_p6) / (1 - csi_S_p4(end));
csi_S_n6 = csi_S_n4(end) * (1 - csi_S_p6) / (1 - csi_S_p4(end));
E6 = E_a + csi_S_p6 * (E_m_p - E_a) + csi_S_n6 * (E_m_n - E_a) + csi_T6 * (E_m_T - E_a);
eps6 = sigma6 ./ E6 + eps_L_p * csi_S_p6 + eps_L_n * csi_S_n6;
sigma_plot = [sigma;sigma2;sigma3;sigma4;sigma5;sigma6];
eps_plot = [eps;eps2;eps3;eps4;eps5;eps6];

%%
figure
plot(eps_plot,sigma_plot)
grid on
