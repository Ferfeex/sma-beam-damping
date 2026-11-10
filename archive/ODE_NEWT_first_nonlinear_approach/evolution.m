function [epsilon,xisp,xism,xiT] = evolution(x0) %x0 vettore spostamento trave
%Material
Ea = 34.9e3;
Emp = 140e3;
Emm = 100e3;
EmT = 70e3;
ssp = 0;
sfp = 120;
ssm = 0;
sfm = -100;
eps_Lp = 0.0012;
eps_Lm = -0.0014;

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
Xa = [x0 ; zeros(2,1)];
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * Xa;
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
    xisp(indp) = 1/2 * (1-cos((sigma_i(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos((sigma_i(indp)-ssp)/(sfp-ssp)));
    xiT(indp) = xiT0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
    xism(indp) = xism0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
    xism(indm) = 1/2 * (1-cos((sigma_i(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos((sigma_i(indm)-ssm)/(sfm-ssm)));
    xiT(indm) = xiT0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
    xisp(indm) = xisp0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
    sigmas = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
    
    F = sigma_i-sigmas;
    index = find(abs(F)>10^-3);
    indp = indp(ismember(indp,index));
    indm = indm(ismember(indm,index));
    
    Dxisp(indp) = 1/2 * pi/(sfp-ssp) * sin(pi*(sigma_i(indp)-ssp)/(sfp-ssp));
    DxiT(indp) = Dxisp(indp) .* (-xiT(indp));
    Dxism(indm) = 1/2 * pi/(sfm-ssm) * sin(pi*(sigma_i(indm)-ssm)/(sfm-ssm));
    DxiT(indm) = Dxism(indm) .* (-xiT(indm));
    
    dE = Dxisp*(Emp-Ea) + Dxism*(Emm-Ea) + DxiT*(EmT-Ea);
    dF = 1 - (dE.*(epsilon-eps_Lp*xisp-eps_Lm*xism) + (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)).*(-eps_Lp*Dxisp-eps_Lm*Dxism));
    sigma_i(index) = sigma_i(index) - F(index)./dF(index);
    i = i+1;
end
xisp=xisp(:);
xism=xism(:);
xiT=xiT(:);
epsilon=epsilon(:);
end
