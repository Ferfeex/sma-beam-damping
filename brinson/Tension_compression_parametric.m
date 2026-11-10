clear all
clc
close all

% Young's
E_a = 34.9e3; %[MPa]
E_m_p = 26e3;%140e3; %[MPa]
E_m_n = 33e3;%90e3;
E_m_T = 26e3;%100e3;
% Other SMA parameters
Teta = 0.55; %[MPa/°C]
Mf = -55; %[°C]
Ms = -35; %[°C]
As = -30; %[°C]
Af = -15; %[°C]
C_M_p = 1; %[MPa/°C]
C_A_p = 3; %[MPa/°C]
C_M_n = 0.2;
C_A_n = 3;
sigma_s_cr_p = 0; %[MPa]
sigma_f_cr_p = 188;%180; %[MPa]
sigma_s_cr_n = 0;
sigma_f_cr_n = -220;%-180;
% [2026 cleanup] thesis Table 3.2 values. The stored 2020 file had 0.0086 / -0.00746, which do not
% reproduce the thesis loss-factor curve (Fig. 2.28).
eps_L_p = 0.016;     % was 0.0086
eps_L_n = -0.0146;   % was -0.00746
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);
T = -60;
% Bekker-Brinson Algorithm
sigma_plot = [];
eps_plot = [];
xisp = 0;
xism = 0;
xiT = 0;
%%%%%%%%%% EPS REVERSAL INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [2026 cleanup] vector of reversal amplitudes 0.1 % ... 2.5 % (thesis Fig. 2.28 abscissa), reconstructed:
% the stored file had the single value 0.025, which makes the polynomial fit degenerate.
eps_upp = 0.001:0.001:0.025;   % was 0.025
%%%%%%%%%% EPS REVERSAL INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_downn = -eps_upp;
Loss_factor = 0;
sigma = 0;
E = E_m_T;
csi_S_p = 0;
csi_S_n = 0;
csi_T = 1;
eps4 = 0;
eu = [];
for kk = eps_upp
%     Inizializzazione
    eps_up = kk;
    eps_down = -eps_up;
    eu = [eu eps_up];
% Cicli all'interno dello stesso ciclo Tensione-Compressione
    for ii = 1:2
% Tension, per eventuale prova con ssp > 0 e/o ssm < 0
%         sigma = linspace(sigma_s_cr_p, sigma_s_cr_p, 1000)';
%         eps = sigma / E(end) + eps_L_p * csi_S_p(end) + eps_L_n * csi_S_n(end);
% Dissipated energy tension
%         E_diss1 = trapz(eps,sigma);
% Reorientation 1
        sigma2 = linspace(0, sigma_f_cr_p, 1000)';
        Y2 = cos(pi * (sigma2 - sigma_s_cr_p) / (sigma_f_cr_p - sigma_s_cr_p));
        csi_S_p2 = (1 - Y2) / 2 + csi_S_p(end) * (1 + Y2) / 2;
        csi_T2 = csi_T(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
        csi_S_n2 = csi_S_n(end) * (1 - csi_S_p2) / (1 - csi_S_p(end));
        E2 = E_a + csi_S_p2 * (E_m_p - E_a) + csi_S_n2 * (E_m_n - E_a) + csi_T2 * (E_m_T - E_a);
        eps2 = sigma2 ./ E2 + eps_L_p * csi_S_p2 + eps_L_n * csi_S_n2;
        eps2 = eps2(eps2 <= eps_up);
        sigma2 = sigma2(1:length(eps2));
        csi_S_p2 = csi_S_p2(1:length(eps2));
        csi_S_n2 = csi_S_n2(1:length(eps2));
        csi_T2 = csi_T2(1:length(eps2));
        E2 = E2(1:length(eps2));
% Energy
        E_diss2 = trapz(eps2,sigma2);
%         U_max2 = 1/2 * eps2(end) * sigma2(end);
% Aggiornamento XI ed E plot
xisp = [xisp ; csi_S_p2];
xism = [xism ; csi_S_n2];
xiT = [xiT ; csi_T2];
E = [E ; E2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if sigma2(end) <= 40
%             sigma_s_cr_n = sigma2(end);
%         else
%             sigma_s_cr_n = 40;
%         end 
% Compression linear
        sigma3 = linspace(sigma2(end), sigma_s_cr_n, 1000)';
        eps3 = sigma3 / E2(end) + eps_L_p * csi_S_p2(end) + eps_L_n * csi_S_n2(end);
% Aggiornamento XI plot
xisp = [xisp ; xisp(end)*ones(length(sigma3),1)];
xism = [xism ; xism(end)*ones(length(sigma3),1)];
xiT = [xiT ; xiT(end)*ones(length(sigma3),1)];
E = [E ; E(end)*ones(length(sigma3),1)];
% Energy
        E_diss3 = trapz(eps3,sigma3);
% Reorientation 2
        sigma4 = linspace(sigma3(end), sigma_f_cr_n, 1000)';
        Y4 = cos(pi * (sigma4 - sigma_s_cr_n) / (sigma_f_cr_n - sigma_s_cr_n));
        csi_S_n4 = (1 - Y4) / 2 + csi_S_n2(end) * (1 + Y4) / 2;
        csi_S_p4 = csi_S_p2(end) * (1 - csi_S_n4) / (1 - csi_S_n2(end));
        csi_T4 = csi_T2(end) * (1 - csi_S_n4) / (1 - csi_S_n2(end));
        E4 = E_a + csi_S_p4 * (E_m_p - E_a) + csi_S_n4 * (E_m_n - E_a) + csi_T4 * (E_m_T - E_a);
        eps4 = sigma4 ./ E4 + eps_L_p * csi_S_p4 + eps_L_n * csi_S_n4;
        eps4 = eps4(eps4 >= eps_down);
        sigma4 = sigma4(1:length(eps4));
        csi_S_p4 = csi_S_p4(1:length(eps4));
        csi_S_n4 = csi_S_n4(1:length(eps4));
        csi_T4 = csi_T4(1:length(eps4));
        E4 = E4(1:length(eps4));
% Energy
        E_diss4 = trapz(eps4,sigma4);
%         U_max4 = 1/2 * eps4(end) .* sigma4(end);
% Aggiornamento XI plot
xisp = [xisp ; csi_S_p4];
xism = [xism ; csi_S_n4];
xiT = [xiT ; csi_T4];
E = [E ; E4];
% Per eventuale prova con sfp variabile con punto di arrivo
%         if sigma4(end) >= -10
%             sigma_s_cr_p = sigma4(end);
%         else
%             sigma_s_cr_p = -10;
%         end        
% Tension linear
        sigma5_lin = linspace(sigma4(end), sigma_s_cr_p, 1000)';
        eps5_lin = sigma5_lin / E4(end) + eps_L_p * csi_S_p4(end) + eps_L_n * csi_S_n4(end);
% Energy 
        E_diss5 = trapz(eps5_lin,sigma5_lin);
% Aggiornamento XI plot
xisp = [xisp ; xisp(end)*ones(length(sigma5_lin),1)];
xism = [xism ; xism(end)*ones(length(sigma5_lin),1)];
xiT = [xiT ; xiT(end)*ones(length(sigma5_lin),1)];
E = [E ; E(end)*ones(length(sigma5_lin),1)];
% Reorientation from ssp < 0 to 0 (particular case)       
%         sigma5 = linspace(sigma5_lin(end), 0, 1000)';
%         Y5 = cos(pi * (sigma5 - sigma_s_cr_p) / (sigma_f_cr_p - sigma_s_cr_p));
%         csi_S_p5 = (1 - Y5) / 2 + csi_S_p4(end) * (1 + Y5) / 2;
%         csi_T5 = csi_T4(end) * (1 - csi_S_p5) / (1 - csi_S_p4(end));
%         csi_S_n5 = csi_S_n4(end) * (1 - csi_S_p5) / (1 - csi_S_p4(end));
%         E5 = E_a + csi_S_p5 * (E_m_p - E_a) + csi_S_n5 * (E_m_n - E_a) + csi_T5 * (E_m_T - E_a);
%         eps5 = sigma5 ./ E5 + eps_L_p * csi_S_p5 + eps_L_n * csi_S_n5;
% Energy
%         E_diss5 = trapz(eps5,sigma5);
% Putting together the variables for the plot
        sigma_plot = [sigma_plot; sigma2; sigma3; sigma4; sigma5_lin];
        eps_plot = [eps_plot; eps2; eps3; eps4; eps5_lin;];        
% Re-inizializzazione variabili
        csi_S_p(end) = csi_S_p4(end);
        csi_S_n(end) = csi_S_n4(end);
        csi_T(end) = csi_T4(end);
% All the energies within the repetition of the same eps-revs cycle
        E_diss_cycle(ii,1) = E_diss2 + E_diss4 + E_diss5;
    end
% Energies of the different Tension-Compression cycles
    E_diss_strain = E_diss_cycle(end);
%     U_max = U_max2 + U_max4;
    U_max = E_diss2 - 1/2*(E_diss2+E_diss3) + E_diss4 - 1/2*(E_diss4+E_diss5);
    Loss_factor(end+1) = E_diss_strain / (2 * pi * U_max)
end
% EPSILON - SIGMA
figure
plot(eps_plot * 100, sigma_plot, 'LineWidth', 1.5);
xlabel('Strain \epsilon [%]')
ylabel('Stress \sigma [MPa]')
grid on
title('Ferro-Elasticity')
% EPSILON_REV - LOSS_FACTOR
p = polyfit(eu,Loss_factor(2:end),4)
% [2026 cleanup] the original never saved p (data/p.mat was saved by hand). It is now written to
% results/p_regen.mat; the SHIPPED data/p.mat (the original one) is NOT overwritten.
% To use the regenerated curve:  cfg = case_base(); cfg.p_file = fullfile('results','p_regen.mat');
here = fileparts(mfilename('fullpath'));
if ~exist(fullfile(here,'..','results'),'dir'), mkdir(fullfile(here,'..','results')); end
save(fullfile(here,'..','results','p_regen.mat'),'p');
xx = linspace(0,max(eps_upp),1000);   % [2026 cleanup] max(): eps_upp is now a vector
yy = p(1)*xx.^4+p(2)*xx.^3+p(3)*xx.^2+p(4)*xx.^1+p(5);
figure
plot(xx*100,yy,'LineWidth',2)
hold on
plot(eu * 100, Loss_factor(2:end),'o')
hold off
grid on
xlabel('Upper Reversal Strain [%]')
ylabel('Loss Factor \eta')
title('Loss Factor Estimation Hysteresis cycles')
% SIGMA - XISP
figure
plot(sigma_plot , xisp(2:end), 'LineWidth', 1.5)
grid on
xlabel('\sigma [MPa]')
ylabel('Tensile Detwinning Fraction \xi_s^+')
title('Evolution of Martensitic tensile detwinning fraction along cycling')
% SIGMA - XISM
figure
plot(sigma_plot , xism(2:end), 'LineWidth', 1.5)
grid on
xlabel('\sigma [MPa]')
ylabel('Compressive Detwinning Fraction \xi_s^-')
title('Evolution of Martensitic compressive detwinning fraction along cycling')
% SIGMA - XIT
figure
plot(sigma_plot , xiT(2:end), 'LineWidth', 1.5)
grid on
xlabel('\sigma [MPa]')
ylabel('Twinned Martensite Fraction \xi_T')
title('Evolution of Twinned Martensite fraction along cycling')
% SIGMA - E
figure
plot(sigma_plot , E(2:end), 'LineWidth', 1.5)
grid on
xlabel('\sigma [MPa]')
ylabel("Young's modulus [MPa]")
title("Evolution of Young's modulus along cycling")
