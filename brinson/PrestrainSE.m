clear all
clc
close all
% CU-AL-ZN-MN properties
Ea = 34.9e3; %[MPa]
Emp = 26e3; %[MPa]
Emm = 33e3;
EmT = 26e3;
Teta = 0.55; %[MPa/°C]
Mf = -55; %[°C]
Ms = -35; %[°C]
As = -30; %[°C]
Af = -15; %[°C]
CMp = 1; %[MPa/°C]
CAp = 3; %[MPa/°C]
CMm = -0.2;
CAm = -3;
ssp = 35; %[MPa]
sfp = 188; %[MPa]
ssm = 0;
sfm = -220;
epsLp = 0.016;
epsLm = -0.0146;
% Phase diagram
T = -60:0.1:50;
Rmfp = CMp * (T - Ms) + sfp; Rmfp = [Rmfp >= sfp] .* Rmfp;
Rmsp = CMp * (T - Ms) + ssp; Rmsp = [Rmsp >= ssp] .* Rmsp;
Rasp = CAp * (T - As); Rasp = [Rasp >= 0] .* Rasp;
Rafp = CAp * (T - Af); Rafp = [Rafp >= 0] .* Rafp;
TT = T(T<=Ms);
% Simulation
T0 = -10;
xisp0 = 0;
xism0 = 0;
xisp = 0;
xism = 0;
sigma0 = 60;
sigma = sigma0;
E = Ea + xisp0*(Emp-Ea) + xism0*(Emm-Ea);
eps = sigma/Ea;
LF = [];
sigma_plot = [];
eps_plot = [];
eps_max = (CMp*(T0-Ms)+ssp) / Emp + epsLp;
eps_min = (CAp*(T0-Af)) / Ea;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TENSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
key = 15;
for k = 0:key
    for ii = 1:1
        % Load linear
        sigma1 = linspace(sigma(end), CMp*(T0-Ms)+ssp ,1000)';
        eps1 = sigma1/E(end) + epsLp*xisp(end);
        % Parameters
        sigma = [sigma ; sigma1];
        eps = [eps ; eps1];
        xisp = [xisp ; zeros(length(sigma1),1)];
        E = [E ; Ea * ones(length(sigma1),1)];
        % DT
        if CMp*(T0-Ms)+sfp-0.9*((CMp*(T0-Ms)+sfp-60)/(key+1))*k > CMp*(T0-Ms)+ssp
            sigma2 = linspace(sigma(end), CMp*(T0-Ms)+sfp-0.9*((CMp*(T0-Ms)+sfp-60)/(key+1))*k ,1000)';
        else
            sigma2 = linspace(sigma(end), CMp*(T0-Ms)+ssp ,1000)';
        end   
        Y = cos(pi * (sigma2-ssp-CMp*(T0-Ms))/(sfp-ssp));
        xisp2 = (1-Y)/2+xisp0*(1+Y)/2;
        E2 = Ea + xisp2*(Emp-Ea);
        eps2 = epsLp*xisp2+sigma2./E2;
        % Parameters
        sigma = [sigma ; sigma2];
        eps = [eps ; eps2];
        E = [E ; E2(end)*ones(length(sigma2),1)];
        xisp = [xisp ; xisp2];
        xisp0 = xisp(end);
        % Energy
        E12 = trapz([eps1;eps2],[sigma1;sigma2]);
        % Unload linear to sigma0
        sigma3 = linspace(sigma(end),sigma0,1000)';
        eps3 = epsLp*xisp(end)+sigma3/E(end);
        % Parameters
        sigma = [sigma ; sigma3];
        eps = [eps ; eps3];
        xisp = [xisp ; xisp(end)*ones(length(sigma3),1)];
        E = [E ; E(end)*ones(length(sigma1),1)];
        % Energy
        E3 = trapz(eps3,sigma3);
        %%%%% Energy Upper %%%%%
        EU = E12 + E3; % E3 dovrebbe essere negativo
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Unload linear to "ssm"
        sigma4 = linspace(sigma(end), CAp*(T0-As) ,1000)';
        eps4 = epsLp*xisp(end)+sigma4/E(end);
        % Parameters
        sigma = [sigma ; sigma4];
        eps = [eps ; eps4];
        xisp = [xisp ; xisp(end)*ones(length(sigma4),1)];
        E = [E ; E(end)*ones(length(sigma4),1)];
        % IT
        if CAp*(T0-Af)+((60-CAp*(T0-Af))/(key+1))*k < CAp*(T0-As)
            sigma5 = linspace(sigma(end), CAp*(T0-Af)+((60-CAp*(T0-Af))/(key+1))*k ,1000)';
        else
            sigma5 = linspace(sigma(end), CAp*(T0-As) ,1000)';
        end
        Y = cos(pi*(T0-As-sigma5/CAp)/(Af-As));
        xisp5 = xisp0/2*(Y+1);
        E5 = Ea + xisp5*(Emp-Ea);
        eps5 = epsLp * xisp5 + sigma5 ./E5;
        % Parameters
        xisp = [xisp ; xisp5];
        E = [E ; E5];
        eps = [eps ; eps5];
        sigma = [sigma ; sigma5];
        xisp0 = xisp(end);
        % Energy
        E45 = trapz([eps4;eps5],[sigma4;sigma5]);
        % Load linear to sigma0
        sigma6 = linspace(sigma(end), sigma0, 1000)';
        eps6 = epsLp*xisp(end) + sigma6/E(end);
        % Parameters
        sigma = [sigma; sigma6];
        eps = [eps ; eps6];
        xisp = [xisp ; xisp(end)];
        E = [E ; E(end)*ones(length(sigma6),1)];
        % Energy
        E6 = trapz(eps6,sigma6);
        %%%%% Energy Lower %%%%%
        EL = E45 + E6;
    end
    sigma_plot(:,k+1) = [sigma1;sigma2;sigma3;sigma4;sigma5;sigma6];
    eps_plot(:,k+1) = [eps1;eps2;eps3;eps4;eps5;eps6];
    E_diss = EU + EL;
    U_max = E12 + E6- 1/2*E_diss;
    LF = [LF E_diss/(2*pi*U_max)];
    % STRESS-STRAIN
    plot(eps_plot(:,k+1),sigma_plot(:,k+1),'Color',[0, 0.4470, 0.7410])
    grid on
    xlabel('\epsilon')
    ylabel('\sigma [MPa]')
    hold on
end
title('Stress-Strain curve')
% LOSS FACTOR
figure
k = 0:key;
ss=CMp*(T0-Ms)+sfp-0.9*((CMp*(T0-Ms)+sfp-60)/(key+1))*k;
p = polyfit(ss,LF,3);
% [2026 cleanup] the original never saved p (data/p_SE.mat was saved by hand). Written to results/p_SE_regen.mat;
% the shipped data/p_SE.mat is NOT overwritten.
here = fileparts(mfilename('fullpath'));
if ~exist(fullfile(here,'..','results'),'dir'), mkdir(fullfile(here,'..','results')); end
save(fullfile(here,'..','results','p_SE_regen.mat'),'p');
Fun = polyval(p,ss);
plot(ss,LF,'o')
hold on
plot(ss,Fun,'Linewidth',1.5)
grid on
ylabel('\eta')
xlabel('Reversal Stress \sigma_{rev} [MPa]')
title('Loss Factor - Reversal Stress')
