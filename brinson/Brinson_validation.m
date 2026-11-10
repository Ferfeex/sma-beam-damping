% Brinson validation with example on the model

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
C_M = 8; %[MPa/°C]
C_A = 13.8; %[MPa/°C]
sigma_s_cr = 100; %[MPa]
sigma_f_cr = 170; %[MPa]
eps_L = 0.067;
a_M = pi / (Ms - Mf);
a_A = pi / (Af - As);

% Initial conditions

eps = 0.005;
csi = 0.046;
E = E_a + csi * (E_m - E_a);
sigma = 128;
T = 20;
T_0 = T + (E * (eps - eps_L * csi) - sigma) / Teta;

T = linspace(20,(C_A * As - Teta * T_0 + E * (eps - eps_L * csi)) / (C_A - Teta), 1000)';
sigma = E * (eps - eps_L * csi) + Teta * (T - T_0);
eps = ones(1000,1) * eps;
R3 = C_A * (T - As);
R3 = [R3 > 0] .* R3;
figure
plot(T,sigma,T,R3, 'r--')
xlabel('T [°C]')
ylabel('Stress [MPa]')
grid on

delta_T = 0.01;

%Tentativo 1

% while sigma(end) > C_A * (T(end) - Af)
%     T(end+1) = T(end) + delta_T;
%     eps(end+1) = eps(end);
%     x = csi(end) / 2 * (cos(a_A * (T(end) - As - sigma(end) / C_A)) + 1);
%     y = @(x) x - csi(end) / 2 * (cos(a_A * (T(end) - As - (sigma(end) + (E_a + x * (E_m - E_a)) * (eps(end) - eps_L * x) + (E_a + csi(end) * (E_m - E_a)) * (eps_L * csi(end) - eps(end)) + Teta * delta_T) / C_A)) + 1);
%     dy = @(x) 1 + csi(end) / 2 * sin(a_A * (T(end) - As - (sigma(end) + (E_a + x * (E_m - E_a)) * (eps(end) - eps_L * x) + (E_a + csi(end) * (E_m - E_a)) * (eps_L * csi(end) - eps(end)) + Teta * delta_T) / C_A)) * a_A / C_A * (eps_L * (sigma(end) + E_a + x * (E_m - E_a)) - (E_m - E_a) * (eps(end) - eps_L * x));
%     [csi(end+1,1),~] = Newton_sempl(x,1000,1e-10,y,dy,2);
%     sigma(end+1) = sigma(end) + (E_a + csi(end) * (E_m - E_a)) * (eps(end) - eps_L * csi(end)) + (E_a + csi(end-1) * (E_m - E_a)) * (eps_L * csi(end-1) - eps(end)) + Teta * delta_T;
% end

% Tentativo 2
% while sigma(end) > C_A * (T(end) - Af)
%     T(end+1) = T(end) + delta_T;
%     xi = csi(end);
%     E = @(xi) E_a + xi * (E_m - E_a);
%     A = @(xi) a_A * (T(end) - As - (E(xi) * (eps(end) - eps_L * xi) + Teta * (T(end) - T_0)) / C_A);
%     dA = @(xi) - a_A / C_A * (-eps_L * E(xi) + (eps(end) - eps_L * xi) * (E_m - E_a));
%     y = @(xi) xi - 1/2 * (cos(A(xi)) + 1);
%     dy = @(xi) 1 - 1/2 * sin(A(xi) * dA(xi));
%     err_norm = 1;
%     while err_norm >= 1e-3
%         xii = xi - 2 * y(xi) / dy(xi);
%         err_norm = norm(xi - xii) / norm(xi);
%         xi = xii;
%     end
%     csi(end+1,1) = xi;
%     sigma(end+1,1) = E(csi(end)) * (eps - eps_L * csi(end)) + Teta * (T(end) - T_0);
% end    

%Tentativo 3 WINNER CAZZO
while sigma(end) > C_A * (T(end) - Af)
    T(end+1,1) = T(end) + delta_T;
    s = sigma(end);

    
    xi = @(s) 0.046 / 2 * (cos(a_A * (T(end) - As - s / C_A)) + 1);
    d_xi = @(s) a_A * 0.046 / (2 * C_A) * sin(a_A * (T(end) - As - s / C_A));
    y = @(s) s - xi(s)^2 * eps_L * (E_a - E_m) - xi(s) * (eps(end) * E_m - E_a * (eps_L + eps(end))) - E_a * eps(end) - Teta * (T(end) - T_0);
    dy = @(s) 1 - d_xi(s) * (2 * eps_L * xi(s) * (E_a - E_m) + E_m * eps(end) - E_a * (eps_L + eps(end)));
    

    err_norm = 1;
    while err_norm >= 1e-6
        ss = s - y(s) / dy(s);
        err_norm = norm(s - ss) / norm(s);
        s = ss;
    end
    
sigma(end+1,1) = s;
csi(end+1,1) = xi(sigma(end));
end

% sigma_0 = 128;
% csi_0 = 0.046;
% eps_0 = 0.005;
% E_0 = E_a + csi_0 * (E_m - E_a);
% K_0 = sigma_0 - E_0 * (eps_0 - csi_0 * eps_L);
% T_0 = 20;
% 
% while sigma(end) > C_A * (T(end) - Af)
%     s = sigma(end);
%     xi = csi(end)

R1 = C_M * (T - Ms) + sigma_f_cr;
R1 = [R1 > sigma_f_cr] .* R1;
R2 = C_M * (T - Ms) + sigma_s_cr;
R2 = [R2 > sigma_s_cr] .* R2;
R3 = C_A * (T - As);
R3 = [R3 > 0] .* R3;
R4 = C_A * (T - Af);
R4 = [R4 > 0] .* R4;
figure
plot(T,sigma, T,R1,'r--', T,R2,'r--', T,R3,'r--', T,R4,'r--')
title('Increase of stress due to constraint on deformation')
xlabel('T [°C]')
ylabel('Stress [MPa]')
legend('stress behavior','phase transition lines')
grid on


% Back to 20°C
%3
T2 = linspace(T(end), (E_a * eps(end) - Teta * T_0 + C_M * Ms - sigma_s_cr) / (C_M - Teta), 1000)';
sigma2 = E_a * eps(end) + Teta * (T2 - T_0);
sigma = [sigma ; sigma2];
T = [T ; T2];
%4
while T(end) > 20
    T(end+1,1) = T(end) - delta_T;
    s = sigma(end);

    
    xi = @(s) 1/2 * cos(pi / (sigma_s_cr - sigma_f_cr) * (s - sigma_f_cr - C_M * (T(end) - Ms))) + 1/2;
    d_xi = @(s) -1/2 * sin(pi / (sigma_s_cr - sigma_f_cr) * (s - sigma_f_cr - C_M * (T(end) - Ms))) * pi / (sigma_s_cr - sigma_f_cr);
    y = @(s) s - xi(s)^2 * eps_L * (E_a - E_m) - xi(s) * (eps(end) * E_m - E_a * (eps_L + eps(end))) - E_a * eps(end) - Teta * (T(end) - T_0);
    dy = @(s) 1 - d_xi(s) * (2 * eps_L * xi(s) * (E_a - E_m) + E_m * eps(end) - E_a * (eps_L + eps(end)));
    

    err_norm = 1;
    while err_norm >= 1e-6
        ss = s - y(s) / dy(s);
        err_norm = norm(s - ss) / norm(s);
        s = ss;
    end
    
sigma(end+1,1) = s;
csi(end+1,1) = xi(sigma(end));
end

R1 = C_M * (T - Ms) + sigma_f_cr;
R1 = [R1 > sigma_f_cr] .* R1;
R2 = C_M * (T - Ms) + sigma_s_cr;
R2 = [R2 > sigma_s_cr] .* R2;
R3 = C_A * (T - As);
R3 = [R3 > 0] .* R3;
R4 = C_A * (T - Af);
R4 = [R4 > 0] .* R4;
figure
plot(T,sigma, T,R1,'r--', T,R2,'r--', T,R3,'r--', T,R4,'r--')
xlabel('T [°C]')
ylabel('Stress [MPa]')
legend('stress behavior','phase transition lines')
grid on
