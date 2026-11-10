% Brinson model (1993)

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


% Initial conditions
%1
T_0 = -55; %[°C]
sigma = linspace(0, sigma_s_cr, 1000)';
eps = sigma / E_m;
eps_s_cr = eps(end);


%2
sigma2 = linspace(sigma_s_cr, sigma_f_cr, 1000)';
csi_S = 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (sigma2 - sigma_f_cr)) + 1/2;
csi_T = 1 - csi_S;
eps2 = eps_L * csi_S + sigma2 / E_m;

eps_rev = eps2(end);

% Loss factor Total reorientation
DW = E_m * eps_s_cr^2 / 2 + 1/E_m * (sigma_f_cr^2 - sigma_s_cr^2) / 2 + eps_L * pi / (2 * (sigma_f_cr - sigma_s_cr)) * 2 * (sigma_f_cr - sigma_s_cr) / pi + E_m * (eps_L^2 - eps_rev^2) / 2 - E_m * eps_L * (eps_L - eps_rev);
W = E_m * eps_s_cr^2 / 2 + 1/E_m * (sigma_f_cr^2 - sigma_s_cr^2) / 2 + eps_L * pi / (2 * (sigma_f_cr - sigma_s_cr)) * 2 * (sigma_f_cr - sigma_s_cr) / pi;
eta = DW / (pi * (W - 1/2 * DW))
eta2 = DW / (2 * pi * W)


figure
plot(sigma2,csi_S,sigma2,csi_T)
legend('\xi_S','\xi_T')
xlabel('Stress [MPa]')
title('Martensitic fractions, Detwinning')
grid on

%3 
sigma3 = linspace(sigma_f_cr, 0, 1000)';
eps3 = eps_L + sigma3 / E_m;

sigma = [sigma; sigma2; sigma3];
eps = [eps; eps2; eps3];

figure
plot(eps,sigma)
grid on
xlabel('Strain')
ylabel('Stress [MPa]')
title('Stress-strain detwinning')

%% STRESS RECOVERY 4 AND 5


% T = linspace(T_0, (E_m * (eps(end) - eps_L * csi_S(end)) + C_A * As - Teta * T_0) / (C_A - Teta) , 1000)';
% sigma4 = E_m * (eps(end) - eps_L * csi_S(end)) + Teta * (T - T_0);
% T = [ones(3000,1) * T_0; T];
% sigma = [sigma ; sigma4];
% 
% E = E_m;
% csi = 1;
% 
% 
% 
% delta_T = 0.05;
% 
% while sigma(end) > C_A * (T(end) - Af) && sigma(end) < 2000
%     T(end+1,1) = T(end) + delta_T;
%     s = sigma(end);
% 
%     
%     xi = @(s) 1 / 2 * (cos(a_A * (T(end) - As - s / C_A)) + 1);
%     d_xi = @(s) a_A / (2 * C_A) * sin(a_A * (T(end) - As - s / C_A));
%     y = @(s) s - xi(s)^2 * eps_L * (E_a - E_m) - xi(s) * (eps(end) * E_m - E_a * (eps_L + eps(end))) - E_a * eps(end) - Teta * (T(end) - T_0);
%     dy = @(s) 1 - d_xi(s) * (2 * eps_L * xi(s) * (E_a - E_m) + E_m * eps(end) - E_a * (eps_L + eps(end)));
%     
% 
%     err_norm = 1;
%     while err_norm >= 1e-6
%         ss = s - y(s) / dy(s);
%         err_norm = norm(s - ss) / norm(s);
%         s = ss;
%     end
%     
% sigma(end+1,1) = s;
% csi(end+1,1) = xi(sigma(end));
% end
% 
% R1 = C_M * (T - Ms) + sigma_f_cr;
% R1 = [R1 > sigma_f_cr] .* R1;
% R2 = C_M * (T - Ms) + sigma_s_cr;
% R2 = [R2 > sigma_s_cr] .* R2;
% R3 = C_A * (T - As);
% R3 = [R3 > 0] .* R3;
% R4 = C_A * (T - Af);
% R4 = [R4 > 0] .* R4;
% TT = -65:0.1:Ms;
% figure
% % plot(T,sigma,'k','LineWidth',1.5)
% % hold on
% plot(T,R1,'r--', T,R2,'r--', T,R3,'r--', T,R4,'r--', TT, ones(length(TT), 1) * sigma_s_cr, 'r--', TT, ones(length(TT), 1) * sigma_f_cr, 'r--')
% xlabel('T [°C]')
% ylabel('Stress [MPa]')
% title('Thermomechanical loading path into phase diagram')
% grid on
% hold off
% 
% ii = length(T) - length(csi) + 1;
% figure
% plot(T(ii:end,1),csi);
% title('Martensitic fraction evolution M --> A')
% xlabel('Temperature [°C]')
% ylabel('\xi')
% grid on

%% STRAIN RECOVERY 


T4 = linspace(T_0, As, 1000)';
eps4 = eps_L - Teta / E_m * (T4 - T_0);



T5 = linspace(As, Af, 1000)';
csi5 = 1/2 * (cos(a_A * (T5 - As)) + 1);
E5 = E_a + csi5 * (E_m - E_a);
eps5 = eps_L * csi5 - Teta / E5 * (T5 - T_0);

T = [ones(3000,1) * T_0 ; T4 ; T5; T_0];
eps = [eps; eps4; eps5; 0];
sigma = [sigma; zeros(2001,1)];


figure
plot3(-eps, T, sigma)
xlabel('\epsilon')
zlabel('\sigma [MPa]')
ylabel('Temperature [°C]')
grid on
title('Shape Memory Effect')

figure
plot(T5,csi5)
title('Total Martensitic fraction, Inverse Transformation')
xlabel('Temperature T [°C]')
ylabel('Total martensitic fraction \xi')
grid on
