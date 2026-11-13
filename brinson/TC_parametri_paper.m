clear all
clc
close all
% Young's
E_a = 34.9e3; %[MPa]
E_m_p = 26e3; %[MPa]
E_m_n = 33e3;
E_m_T = 26e3;
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
sigma_s_cr_p = 35; %[MPa]
sigma_f_cr_p = 180; %[MPa]
sigma_s_cr_n = 0;
sigma_f_cr_n = -220;
eps_L_p = 0.016;
eps_L_n = -0.0146;
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
eps_upp = 0.023;
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
for kk = eps_upp : -0.001 : 0.002
%     Inizializzazione
    eps_up = kk;
    eps_down = -eps_up;
    eu = [eu eps_up];
% Cicli all'interno dello stesso ciclo Tensione-Compressione
    for ii = 1:10
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
% Aggiornamento XI ed E plot
xisp = [xisp ; csi_S_p2];
xism = [xism ; csi_S_n2];
xiT = [xiT ; csi_T2];
E = [E ; E2];
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
% Aggiornamento XI plot
xisp = [xisp ; csi_S_p4];
xism = [xism ; csi_S_n4];
xiT = [xiT ; csi_T4];
E = [E ; E4];
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
    U_max = E_diss2 - 1/2*(E_diss2+E_diss3) + E_diss4 - 1/2*(E_diss4+E_diss5);
    Loss_factor(end+1) = E_diss_strain / (2 * pi * U_max)
end
% EPSILON - SIGMA
figure
plot(eps_plot * 100, sigma_plot);
xlabel('Strain \epsilon [%]')
ylabel('Stress \sigma [MPa]')
grid on
title('Stress-Strain diagram')
% EPSILON_REV - LOSS_FACTOR
figure
plot(eu * 100, Loss_factor(2:end),'o-')
grid on
xlabel('Upper Reversal Strain [%]')
ylabel('Loss Factor \eta')
title('Loss Factor Estimation Hysteresis cycles')
% SIGMA - XISP
figure
plot(sigma_plot , xisp(2:end))
grid on
xlabel('\sigma [MPa]')
ylabel('Tensile Detwinning Fraction \xi_s^+')
title('Evolution of Martensitic tensile detwinning fraction along cycling')
% SIGMA - XISM
figure
plot(sigma_plot , xism(2:end))
grid on
xlabel('\sigma [MPa]')
ylabel('Compressive Detwinning Fraction \xi_s^-')
title('Evolution of Martensitic compressive detwinning fraction along cycling')
% SIGMA - XIT
figure
plot(sigma_plot , xiT(2:end))
grid on
xlabel('\sigma [MPa]')
ylabel('Twinned Martensite Fraction \xi_T')
title('Evolution of Twinned Martensite fraction along cycling')
% SIGMA - E
figure
plot(sigma_plot , E(2:end))
grid on
xlabel('\sigma [MPa]')
ylabel("Young's modulus [MPa]")
title("Evolution of Young's modulus along cycling")
