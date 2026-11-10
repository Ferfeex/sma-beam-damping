clear
clc
close all

%Material
Ea = 34.9e3;
Emp = 140e3;
Emm = 100e3;
EmT = 70e3;
ssp = 0;
sfp = 80;
ssm = 0;
sfm = -70;
eps_Lp = 0.0012;
eps_Lm = -0.0014;

% Initial displacement field (and null velocity field)
% XA = [[1.68572689366932,0.0249698216568307,3.30823750726361,0.0235474570906879,4.80431602126453,0.0211768609379581,6.11074752033265,0.0178580499400841,7.16431840296102,0.0135910466426254,7.90181673991894,0.00837587849656379,8.26003256394894,0.00221257669348144,8.19550982388027,-0.00400999647101798,7.74380173678819,-0.00940297806783367,6.96021551984661,-0.0139663350252603,5.90006051751249,-0.0177000368267950,4.61864798770562,-0.0206040567766679,3.17129081132115,-0.0226783730347844,1.61330314037335,-0.0239229694230469,0.0254439443373185,-0.0243378360029651]';zeros(2,1)];
% X0s = XA(1:30);
X0s=[0.00421475775365677,6.24310731598045e-05,0.00827145710652474,5.88747648352622e-05,0.0120120411683463,5.29476218662368e-05,0.0152784560125262,4.46496991043185e-05,0.0179126520146303,3.39810703338018e-05,0.0197565850193508,2.09418252740778e-05,0.0206522172795097,5.53206574024694e-06,0.0204909004589998,-1.00258933824359e-05,0.0193615229942669,-2.35097351619935e-05,0.0174023630045132,-3.49193516989787e-05,0.0147517055525488,-4.42546433352411e-05,0.0115478419531615,-5.15155228007974e-05,0.00792906882912855,-5.67019186347819e-05,0.00403368696378613,-5.98137778655045e-05,6.36165130682430e-05,-6.08510679362297e-05]';

% Strain field
N = 16; % Numero nodi
L = 1000; % [mm] lunghezza trave
Lk = L/(N-1);
Extr = zeros(4, 2*N, N-1);
for k=2:N-2
    Extr(1:4,2*k-3:2*k,k) = eye(4);
end
% Extraction matrices for boundary elements, hinges on both sides
Extr(3:4,1:2,1) = eye(2);
Extr(2, 2*N - 3, 1) = 1;
Extr(1, 2*N-1, 1) = 1;
Extr(1:2,2*N-5:2*N-4,N-1) = eye(2);
Extr(3,2*N,N-1) = 1;
Extr(4,2*N-2,N-1) = 1;
Extr_tr = permute(Extr, [2,1,3]); % Trasposing every extraction matrix

csi = linspace(0,Lk,5)';
y = linspace(10,-10,10)';
[CSI,Y] = meshgrid(csi,y);
Bk = zeros(length(csi)*length(y),4);
Bk1 = -Y .* (12*CSI/Lk^3-6/Lk^2);
Bk(:,1) = Bk1(:);
Bk2 = -Y .* (6/Lk^2*CSI - 4/Lk);
Bk(:,2) = Bk2(:);
Bk(:,3) = -Bk(:,1);
Bk4 = -Y .* (6/Lk^2*CSI - 2/Lk);
Bk(:,4) = Bk4(:);
Bk = repmat({Bk},1,N-1);
B = blkdiag(Bk{:});
X_k = zeros(4*(N-1),1);
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
end
epsilons = B * X_k;
epsilon = reshape(epsilons,length(y),[]);

sigma_i = EmT .* epsilon;
xisp = zeros(size(sigma_i));
xism = xisp;
xiT = ones(size(sigma_i));

index = 1;
indp = find(sigma_i>0);
indm = find(sigma_i<0);
i = 0;
while ~isempty(index)
    
    xisp(indp) = 1/2 * (1-cos((sigma_i(indp)-ssp)/(sfp-ssp)));
    xiT(indp) = 1 - xisp(indp);
    xism(indm) = 1/2 * (1-cos((sigma_i(indm)-ssm)/(sfm-ssm)));
    xiT(indm) = 1 - xism(indm);
    sigmas = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
    
    F = sigma_i-sigmas;
    index = find(abs(F)>10^-6);
    indp = intersect(index,indp);
    indm = intersect(index,indm);
    
    dxisp = zeros(size(xisp)); dxism = zeros(size(xism)); dxiT = zeros(size(xiT));
    dxisp(indp) = 1/2 * pi/(sfp-ssp) * sin(pi*(sigma_i(indp)-ssp)/(sfp-ssp));
    dxiT(indp) = dxisp(indp) .* (-xiT(indp));
    dxism(indm) = 1/2 * pi/(sfm-ssm) * sin(pi*(sigma_i(indm)-ssm)/(sfm-ssm));
    dxiT(indm) = dxism(indm) .* (-xiT(indm));
    
    dE = dxisp*(Emp-Ea) + dxism*(Emm-Ea) + dxiT*(EmT-Ea);
    dF = 1 - (dE.*(epsilon-eps_Lp*xisp-eps_Lm*xism) + (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)).*(-eps_Lp*dxisp-eps_Lm*dxism));
    sigma_i(index) = sigma_i(index) - F(index)./dF(index);
    i = i+1;
end
i

x = linspace(0,L,75)';
figure
mesh(x,y,sigma_i)
figure
subplot(3,1,1)
mesh(x,y,xisp)
subplot(3,1,2)
mesh(x,y,xism)
subplot(3,1,3)
mesh(x,y,xiT)

IN_COND = [X0s; xisp(:); xism(:); xiT(:); sigma_i(:); epsilons];
