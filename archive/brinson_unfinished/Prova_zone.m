clear
clc
close all

% Material properties
E_a = 67e3; %[MPa]
E_m = 26.3e3; %[MPa]
Teta = 0.55; %[MPa/°C]
Mf = 9; %[°C]
Ms = 18.4; %[°C]
As = 34.5; %[°C]
Af = 49; %[°C]
C_M = 10.8; %[MPa/°C]
C_A = 13.8; %[MPa/°C]
sigma_s_cr = 100; %[MPa]
sigma_f_cr = 170; %[MPa]
eps_L = 0.067;
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);

sigma = [50, 137, 200];
Temp = [5, 13, 20, 25, 30, 50, 75];

% Trials to verify zone function works

zero = Zone( sigma(1), Temp(1), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
t = Zone( sigma(1), Temp(2), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(1), Temp(3), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(1), Temp(5), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
A = Zone( sigma(1), Temp(6), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(1), Temp(7), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
o = Zone( sigma(2), Temp(1), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
ot = Zone( sigma(2), Temp(2), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
M = Zone( sigma(2), Temp(3), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(2), Temp(5), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
A = Zone( sigma(2), Temp(6), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(2), Temp(7), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(3), Temp(1), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
M = Zone( sigma(3), Temp(4), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(3), Temp(5), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
A = Zone( sigma(3), Temp(6), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)
zero = Zone( sigma(3), Temp(7), Mf, Ms, As, Af, sigma_s_cr, sigma_f_cr, C_M, C_A)


T = linspace(0, 80, 1000);
R1 = C_M * (T - Ms) + sigma_f_cr;
R1 = [R1 > sigma_f_cr] .* R1;
R2 = C_M * (T - Ms) + sigma_s_cr;
R2 = [R2 > sigma_s_cr] .* R2;
R3 = C_A * (T - As);
R3 = [R3 > 0] .* R3;
R4 = C_A * (T - Af);
R4 = [R4 > 0] .* R4;
TT = [T <= Ms] .* T;
figure
plot(T,R1,'r--', T,R2,'r--', T,R3,'r--', T,R4,'r--', TT, ones(length(TT), 1) * sigma_s_cr, 'r--', TT, ones(length(TT), 1) * sigma_f_cr, 'r--')
xlabel('T [°C]')
ylabel('Stress [MPa]')
grid on
