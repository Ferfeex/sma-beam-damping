function y0 = inizcon(x0) %x0 vettore spostamento trave
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

y0 = [zeros(30,1); x0; xisp(:); xism(:); xiT(:); sigma_i(:); epsilons];
end
