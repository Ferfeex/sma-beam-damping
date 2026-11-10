function [epsilon,sigma,xisp,xism,xiT] = evolution(x0,xisp0,xism0,xiT0,xisp_in,xism_in,xiT_in,sigma_old) %x0 vettore spostamento trave
%Material
Ea = 34.9e3;
Emp = 30e3;
Emm = 30e3;
EmT = 30e3;
ssp = 0;
sfp = 188;
ssm = 0;
sfm = -220;
eps_Lp = 0.012;
eps_Lm = -0.0146;

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
Bk = zeros(length(csi)*length(y),4);
[CSI,Y] = meshgrid(csi,y);
Bk1 = -Y .* (12*CSI/Lk^3-6/Lk^2);
Bk(:,1) = Bk1(:);
Bk2 = -Y .* (6*CSI/Lk^2 - 4/Lk);
Bk(:,2) = Bk2(:);
Bk(:,3) = -Bk1(:);
Bk4 = -Y .* (6*CSI/Lk^2 - 2/Lk);
Bk(:,4) = Bk4(:);
Bk = repmat({Bk},1,N-1);
B = blkdiag(Bk{:});
X_k = zeros(4*(N-1),1);
Xa = [x0' ; zeros(2,1)];
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * Xa;
end
epsilon = B * X_k;
E = Ea + xisp_in*(Emp-Ea) + xism_in*(Emm-Ea) + xiT_in*(EmT-Ea);
sigma0 = E .* epsilon;
X0 = [sigma0 ; xisp_in ; xism_in ; xiT_in];
JB = zeros(3000);
EY = repmat(eye(750),3,3);
EY2 = repmat(eye(750),1,3);
JB(1:2250,1:2250) = EY;
JB(2251:3000,751:3000) = EY2;
JB(1:750,2251:3000) = eye(750);
options = optimoptions('fsolve','JacobPattern',JB,'Algorithm','trust-region');
X = fsolve(@(X)myfunsolve(X,epsilon,xisp0,xism0,xiT0,sigma_old,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm) , X0 , options);
sigma = X(1:750);
xisp = X(751:1500);
xism = X(1501:2250);
xiT = X(2251:3000);
end
