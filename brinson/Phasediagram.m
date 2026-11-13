clear
clc
close all

% Material properties
E_a = 34.9e3; %[MPa]
E_m = 26.3e3; %[MPa]
Teta = 0.55; %[MPa/°C]
Mf = -55; %[°C]
Ms = -35; %[°C]
As = -30; %[°C]
Af = 15; %[°C]
C_M = 1; %[MPa/°C]
C_A = 3; %[MPa/°C]
sigma_s_cr = 40; %[MPa]
sigma_f_cr = 188; %[MPa]
eps_L = 0.016;
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);

T = -55:0.01:50;

R1 = C_M * (T - Ms) + sigma_f_cr;
R1 = [R1 > sigma_f_cr] .* R1;
R2 = C_M * (T - Ms) + sigma_s_cr;
R2 = [R2 > sigma_s_cr] .* R2;
R3 = C_A * (T - As);
R3 = [R3 > 0] .* R3;
R4 = C_A * (T - Af);
R4 = [R4 > 0] .* R4;

figure
plot(T,R1,T,R2,T,R3,T,R4)
grid on
xlabel('Temperature [°C]')
ylabel('Stress \sigma [MPa]')

T = 28.5;

sigma = linspace(0, C_M * (T - Ms) + sigma_s_cr, 1000)';
eps = sigma / E_a;
csi = zeros(1000,1);

sigma2 = linspace(sigma(end), C_M * (T - Ms) + sigma_f_cr, 1000)';
csi2 = 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma2 - sigma_f_cr - C_M * (T - Ms))) + 1/2;
E2 = E_a + csi2 * (E_m - E_a);
eps2 = eps_L * csi2 + sigma2 ./ E2;

sigma3 = linspace(sigma2(end), C_A * (T - As), 1000)';
eps3 = eps_L + sigma3 / E2(end);
csi3 = ones(1000,1);

sigma4 = linspace(sigma3(end), C_A * (T - Af), 1000)';
csi4 = 1/2 * (cos(a_A * (T - As - sigma4 / C_A)) + 1);
E4 = E_a + csi4 * (E_m - E_a);
eps4 = eps_L * csi4 + sigma4 ./E4;

sigma5 = linspace(sigma4(end), 0, 1000)';
eps5 = sigma5 / E_a;
csi5 = zeros(1000,1);

sigma6 = linspace(sigma(end), 190.5, 1000)';
csi6 = 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma6 - sigma_f_cr - C_M * (T - Ms))) + 1/2;
E6 = E_a + csi6 * (E_m - E_a);
eps6 = eps_L * csi6 + sigma6 ./ E6;

sigma7 = linspace(sigma6(end), C_A * (T - As), 1000)';
eps7 = eps_L * csi6(end) + sigma7 / E6(end);
csi7 = ones(1000,1) * csi6(end);

sigma8 = linspace(sigma7(end), 110, 1000)';
csi8 = csi7(end)/2 * (cos(a_A * (T - As - sigma8 / C_A)) + 1);
E8 = E_a + csi8 * (E_m - E_a);
eps8 = eps_L * csi8 + sigma8 ./E8;

% At 110 MPa the point is already inside the martensite transformation zone (110 > C_M (T - Ms) + sigma_s_cr).

sigma9 = linspace(sigma8(end), C_M * (T - Ms) + sigma_f_cr, 1000)';
csi9 = (1 - csi8(end)) / 2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma9 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi8(end)) / 2;
E9 = E_a + csi9 * (E_m - E_a);
eps9 = eps_L * csi9 + sigma9 ./ E9;

sigma10 = linspace(sigma9(end), C_A * (T - As), 1000)';
eps10 = eps_L * csi9(end) + sigma10 / E9(end);
csi10 = ones(1000,1) * csi9(end);

sigma11 = linspace(sigma10(end), 110, 1000)';
csi11 = csi10(end)/2 * (cos(a_A * (T - As - sigma11 / C_A)) + 1);
E11 = E_a + csi11 * (E_m - E_a);
eps11 = eps_L * csi11 + sigma11 ./E11;

sigma12 = linspace(sigma11(end), sigma6(end), 1000)';
csi12 = (1 - csi11(end)) / 2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma12 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi11(end)) / 2;
E12 = E_a + csi12 * (E_m - E_a);
eps12 = eps_L * csi12 + sigma12 ./ E12;

sigma13 = linspace(sigma9(end), C_A * (T - As), 1000)';
eps13 = eps_L * csi12(end) + sigma13 / E12(end);
csi13 = ones(1000,1) * csi12(end);

sigma14 = linspace(sigma13(end), 110, 1000)';
csi14 = csi13(end)/2 * (cos(a_A * (T - As - sigma14 / C_A)) + 1);
E14 = E_a + csi14 * (E_m - E_a);
eps14 = eps_L * csi14 + sigma14 ./E14;

sigma15 = linspace(sigma14(end), sigma6(end), 1000)';
csi15 = (1 - csi14(end)) / 2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma15 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi14(end)) / 2;
E15 = E_a + csi15 * (E_m - E_a);
eps15 = eps_L * csi15 + sigma15 ./ E15;

sigma16 = linspace(sigma9(end), C_A * (T - As), 1000)';
eps16 = eps_L * csi15(end) + sigma16 / E15(end);
csi16 = ones(1000,1) * csi15(end);

sigma17 = linspace(sigma16(end), 110, 1000)';
csi17 = csi16(end)/2 * (cos(a_A * (T - As - sigma17 / C_A)) + 1);
E17 = E_a + csi17 * (E_m - E_a);
eps17 = eps_L * csi17 + sigma17 ./E17;

sigma18 = linspace(sigma17(end), sigma6(end), 1000)';
csi18 = (1 - csi17(end)) / 2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma18 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi17(end)) / 2;
E18 = E_a + csi18 * (E_m - E_a);
eps18 = eps_L * csi18 + sigma18 ./ E18;

sigma19 = linspace(sigma9(end), C_A * (T - As), 1000)';
eps19 = eps_L * csi18(end) + sigma19 / E18(end);
csi19 = ones(1000,1) * csi18(end);

sigma20 = linspace(sigma19(end), 110, 1000)';
csi20 = csi19(end)/2 * (cos(a_A * (T - As - sigma20 / C_A)) + 1);
E20 = E_a + csi20 * (E_m - E_a);
eps20 = eps_L * csi20 + sigma20 ./E20;

sigma21 = linspace(sigma20(end), sigma6(end), 1000)';
csi21 = (1 - csi20(end)) / 2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma21 - sigma_f_cr - C_M * (T - Ms))) + (1 + csi20(end)) / 2;
E21 = E_a + csi21 * (E_m - E_a);
eps21 = eps_L * csi21 + sigma21 ./ E21;

sigma = [sigma; sigma2; sigma3; sigma4; sigma5; sigma6; sigma7; sigma8; sigma9; sigma10; sigma11; sigma12;sigma13;sigma14;sigma15;sigma16;sigma17;sigma18;sigma19;sigma20;sigma21];
eps = [eps; eps2; eps3; eps4; eps5; eps6; eps7; eps8; eps9; eps10; eps11; eps12;eps13;eps14;eps15;eps16;eps17;eps18;eps19;eps20;eps21];
 
figure
plot(eps, sigma)
grid on
xlabel('\epsilon')
ylabel('\sigma [MPa]')
