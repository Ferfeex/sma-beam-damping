clear all
clc
close all
global k LF
LF = [];
for k = 0
clear
clc
global k LF
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
Rmfm = CMm * (T - Ms) + sfm; Rmfm = [Rmfm <= sfm] .* Rmfm;
Rmsm = CMm * (T - Ms) + ssm; Rmsm = [Rmsm <= ssm] .* Rmsm;
Rasm = CAm * (T - As); Rasm = [Rasm <= 0] .* Rasm;
Rafm = CAm * (T - Af); Rafm = [Rafm <= 0] .* Rafm;
TT = T(T<=Ms);
% Simulation
T0 = -10;
xisp0 = 0;
xism0 = 0;
xisp = 0;
xism = 0;
sigma = 0;
eps = 0;
E = Ea + xisp0*(Emp-Ea) + xism0*(Emm-Ea);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TENSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load linear
sigma2 = linspace(0, CMp * (T0 - Ms) + ssp, 1000)';
eps2 = sigma2 / E(end);
% XI ed E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma2),1)];
xism = [xism ; xism(end)*ones(length(sigma2),1)];
E = [E ; E(end)*ones(length(sigma2),1)];
sigma = [sigma ; sigma2];
eps = [eps ; eps2];
% Energy
Et2 = trapz(eps2,sigma2);
% Transformation direct
sigma3 = linspace(sigma(end), CMp * (T0 - Ms) + sfp - 10*k, 1000)';%%%%%%%%%%
Y = cos(pi * (sigma3-ssp-CMp*(T0-Ms))/(sfp-ssp));
xisp3 = (1-Y)/2+xisp0*(1+Y)/2;
E3 = Ea + xisp3*(Emp-Ea);
eps3 = epsLp * xisp3 + sigma3 ./ E3;
% XI E eps sigma plot
xisp = [xisp ; xisp3];
xism = [xism ; xism(end)*ones(length(sigma3),1)];
E = [E ; E3];
eps = [eps ; eps3];
sigma = [sigma ; sigma3];
% AGGIORNAMENTO XI0
xisp0 = xisp(end);
xism0 = xism(end);
xi0 = xisp0 + xism0;
% Energy
Et3 = trapz(eps3,sigma3);
Et23 = Et2 + Et3;
% Unload linear 1
sigma4 = linspace(sigma(end), CAp * (T0 - As), 1000)';
eps4 = epsLp*xisp(end) + sigma4 / E(end);
% XI E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma4),1)];
xism = [xism ; xism(end)*ones(length(sigma4),1)];
E = [E ; E(end)*ones(length(sigma4),1)];
eps = [eps ; eps4];
sigma = [sigma ; sigma4];
% Transformation inverse
sigma5 = linspace(sigma(end), CAp * (T0 - Af), 1000)';
Y = cos(pi*(T0-As-sigma5/CAp)/(Af-As));
xi = xi0/2*(Y+1);
xisp5 = xisp0/xi0*xi;
E5 = Ea + xisp5*(Emp-Ea);
eps5 = epsLp * xisp5 + sigma5 ./E5;
% XI E eps sigma plot
xisp = [xisp ; xisp5];
xism = [xism ; xism(end)*ones(length(sigma5),1)];
E = [E ; E5];
eps = [eps ; eps5];
sigma = [sigma ; sigma5];
% AGGIORNAMENTO XI0
xisp0 = xisp(end);
xism0 = xism(end);
xi0 = xisp0 + xism0;
% Unload linear 2
sigma6 = linspace(sigma(end), 0, 1000)';
eps6 = sigma6 / E(end);
% XI E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma6),1)];
xism = [xism ; xism(end)*ones(length(sigma6),1)];
E = [E ; E(end)*ones(length(sigma6),1)];
eps = [eps ; eps6];
sigma = [sigma ; sigma6];
%Energy
Et = trapz(eps,sigma);
Umaxt = Et23 - 1/2*Et;
LFt = Et / (2*pi*Umaxt);
%%%%%%%%%%%%%%%%%%%%%%%%% COMPRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load linear
sigma2 = linspace(0, CMm * (T0 - Ms) + ssm, 1000)';
eps2 = sigma2 / E(end);
% XI ed E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma2),1)];
xism = [xism ; xism(end)*ones(length(sigma2),1)];
E = [E ; E(end)*ones(length(sigma2),1)];
sigma = [sigma ; sigma2];
eps = [eps ; eps2];
% Energy
Ec2 = trapz(eps2,sigma2);
% Transformation direct
sigma3 = linspace(sigma(end), CMm * (T0 - Ms) + sfm + 10*k, 1000)';%%%%%%%%%%
Y = cos(pi * (sigma3-ssm-CMm*(T0-Ms))/(sfm-ssm));
xism3 = (1-Y)/2+xism0*(1+Y)/2;
E3 = Ea + xism3*(Emm-Ea);
eps3 = epsLm * xism3 + sigma3 ./ E3;
% XI E eps sigma plot
xism = [xism ; xism3];
xisp = [xisp ; xisp(end)*ones(length(sigma3),1)];
E = [E ; E3];
eps = [eps ; eps3];
sigma = [sigma ; sigma3];
% AGGIORNAMENTO XI0
xisp0 = xisp(end);
xism0 = xism(end);
xi0 = xisp0 + xism0;
% Energy
Ec3 = trapz(eps3,sigma3);
% Unload linear 1
sigma4 = linspace(sigma(end), CAm * (T0 - As), 1000)';
eps4 = epsLm*xism(end) + sigma4 / E(end);
% XI E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma4),1)];
xism = [xism ; xism(end)*ones(length(sigma4),1)];
E = [E ; E(end)*ones(length(sigma4),1)];
eps = [eps ; eps4];
sigma = [sigma ; sigma4];
% Energy
Ec4 = trapz(eps4,sigma4);
% Transformation inverse
sigma5 = linspace(sigma(end), CAm * (T0 - Af), 1000)';
Y = cos(pi*(T0-As-sigma5/CAm)/(Af-As));
xi = xi0/2*(Y+1);
xism5 = xism0/xi0*xi;
E5 = Ea + xism5*(Emm-Ea);
eps5 = epsLm * xism5 + sigma5 ./E5;
% XI E eps sigma plot
xism = [xism ; xism5];
xisp = [xisp ; xisp(end)*ones(length(sigma5),1)];
E = [E ; E5];
eps = [eps ; eps5];
sigma = [sigma ; sigma5];
% AGGIORNAMENTO XI0
xisp0 = xisp(end);
xism0 = xism(end);
xi0 = xisp0 + xism0;
% Energy
Ec5 = trapz(eps5,sigma5);
% Unload linear 2
sigma6 = linspace(sigma(end), 0, 1000)';
eps6 = sigma6 / E(end);
% XI E eps sigma plot
xisp = [xisp ; xisp(end)*ones(length(sigma6),1)];
xism = [xism ; xism(end)*ones(length(sigma6),1)];
E = [E ; E(end)*ones(length(sigma6),1)];
eps = [eps ; eps6];
sigma = [sigma ; sigma6];
% Energy
Ec6 = trapz(eps6,sigma6);
Ec = Ec2 + Ec3 + Ec4 + Ec5 + Ec6;
Umaxc = (Ec2 + Ec3) - 1/2*Ec;
LFc = Ec / (2*pi*Umaxc);
LF = [LF LFt+LFc];
% SIGMA - XISP
figure(1)
subplot(1,2,1)
plot(sigma , xisp,'Color',[0, 0.4470, 0.7410])
grid on
xlabel('\sigma [MPa]')
ylabel('Tensile Detwinning Fraction \xi_s^+')
title('Evolution of Martensitic tensile detwinning fraction')
hold on
% SIGMA - XISM
subplot(1,2,2)
plot(sigma , xism ,'Color',[0, 0.4470, 0.7410])
grid on
xlabel('\sigma [MPa]')
ylabel('Compressive Detwinning Fraction \xi_s^-')
title('Evolution of Martensitic compressive detwinning fraction')
hold on
% SIGMA - E
figure(2)
plot(sigma , E,'Color',[0, 0.4470, 0.7410])
grid on
xlabel('\sigma [MPa]')
ylabel("Young's modulus [MPa]")
title("Evolution of Young's modulus")
hold on
% TEMPERATURE - SIGMA
figure(3)
plot(T,Rmfp,'b',T,Rmsp,'b',T,Rasp,'r',T,Rafp,'r',T,Rmfm,'b',T,Rmsm,'b',T,Rasm,'r',T,Rafm,'r',TT,ssp*ones(length(TT)),'b',TT,sfp*ones(length(TT)),'b',TT,ssm*ones(length(TT)),'b',TT,sfm*ones(length(TT)),'b')
grid on
xlabel('Temperature [°C]')
ylabel('Stress \sigma [MPa]')
title('SMA Phase Diagram')
hold on
% EPSILON - SIGMA
figure(4)
plot(eps * 100, sigma,'Color',[0, 0.4470, 0.7410]);
xlabel('Strain \epsilon [%]')
ylabel('Stress \sigma [MPa]')
title('Super-Elasticity')
hold on
end
%SIGMA_REV-LF
Sigma = [];
for k = 0:12
Sigma = [Sigma 208-10*k -224+10*k];
end
figure
plot(1:length(Sigma),Sigma)
ylabel('\sigma [MPa]')
xlabel('Cycles')
title('Loading cycles')
grid on
figure
plot(Sigma(1:2:end),LF,'o-')
grid on
xlabel('\sigma_{max} [MPa]')
ylabel('Loss Factor \eta')
title('Loss Factor variation with cycles')
