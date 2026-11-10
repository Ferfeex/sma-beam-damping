function bm = build_beam(cfg, E_K)
% BUILD_BEAM  Euler-Bernoulli beam, simply supported, N nodes (N-1 elements).
%   Mass matrix, linear stiffness (modulus E_K), extraction matrices, and the
%   constant B and BB.. matrices used by the nonlinear stiffness integration.
%   This is the geometry/FEM block that was copy-pasted at the top of every
%   FEM1D.m in the original code (variable names unchanged).
N = cfg.geom.N;
L = cfg.geom.L;          % Length [mm]
b = cfg.geom.b;          % Width [mm]
thick = cfg.geom.thick;  % Thickness [mm]
rho = cfg.mat.rho;       % Density [ton/mm^3]
Area = b * thick;
J = 1/12 * b * thick^3;  % Inertia
Lk = L/(N-1);
csi = linspace(0,Lk,cfg.geom.n_csi)';           % x discretization within single element
y = linspace(thick/2,-thick/2,cfg.geom.n_y)';   % y discretization (through the thickness)
% ELEMENT MASS
Mk = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
E = E_K;
Kk = E*J * [  +12./Lk^3  6./Lk^2     -12./Lk^3  6./Lk^2
               6./Lk^2  +4./Lk       -6./Lk^2  +2./Lk
               -12./Lk^3  -6./Lk^2   +12./Lk^3  -6./Lk^2
               6./Lk^2  +2./Lk       -6./Lk^2  +4./Lk ] ;
% EXTRACTION MATRIX
Extr = zeros(4, 2*N, N-1);
for k=2:N-2
    Extr(1:4,2*k-3:2*k,k) = eye(4);
end
Extr(3:4,1:2,1) = eye(2);
Extr(2, 2*N - 3, 1) = 1;
Extr(1, 2*N-1, 1) = 1;
Extr(1:2,2*N-5:2*N-4,N-1) = eye(2);
Extr(3,2*N,N-1) = 1;
Extr(4,2*N-2,N-1) = 1;
Extr_tr = permute(Extr, [2,1,3]);
% ASSEMBLY
M = zeros(2*N,2*N,N-1);
K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
    K(:,:,k) = Extr_tr(:,:,k) * Kk * Extr(:,:,k);
end
M = sum(M,3);
K = sum(K,3);
MFF = M(1:2*N-2,1:2*N-2);
KFF = K(1:2*N-2,1:2*N-2);
% CALCULATION OF B AND B' CONSTANT MATRICES FOR K NONLINEAR MATRIX
%%%%%%%%%%%% 1 %%%%%%%%%%%%
[CSI,Y] = meshgrid(csi,y);
Bk = zeros(length(csi)*length(y),4);
Bk1 = -Y .* (12*CSI/Lk^3-6/Lk^2);
Bk(:,1) = Bk1(:);
Bk2 = -Y .* (6*CSI/Lk^2 - 4/Lk);
Bk(:,2) = Bk2(:);
Bk(:,3) = -Bk1(:);
Bk4 = -Y .* (6*CSI/Lk^2 - 2/Lk);
Bk(:,4) = Bk4(:);
Bk = repmat({Bk},1,N-1);
B = blkdiag(Bk{:});
%%%%%%%%%%%% 2 %%%%%%%%%%%%
BB11 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).^2;
BB12 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).*(-Y .* (6*CSI/Lk^2 - 4/Lk));
BB13 = -BB11;
BB14 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).*(-Y .* (6*CSI/Lk^2 - 2/Lk));
BB21 = BB12;
BB22 = (-Y .* (6*CSI/Lk^2 - 4/Lk)).^2;
BB23 = -BB12;
BB24 = (-Y .* (6*CSI/Lk^2 - 4/Lk)).*(-Y .* (6*CSI/Lk^2 - 2/Lk));
BB31 = BB13;
BB32 = BB23;
BB33 = BB11;
BB34 = -BB14;
BB41 = BB14;
BB42 = BB24;
BB43 = BB34;
BB44 = (-Y .* (6*CSI/Lk^2 - 2/Lk)).^2;

bm = struct('N',N,'L',L,'b',b,'thick',thick,'Lk',Lk,'csi',csi,'y',y, ...
    'Extr',Extr,'Extr_tr',Extr_tr,'MFF',MFF,'KFF',KFF,'B',B, ...
    'BB11',BB11,'BB12',BB12,'BB13',BB13,'BB14',BB14, ...
    'BB21',BB21,'BB22',BB22,'BB23',BB23,'BB24',BB24, ...
    'BB31',BB31,'BB32',BB32,'BB33',BB33,'BB34',BB34, ...
    'BB41',BB41,'BB42',BB42,'BB43',BB43,'BB44',BB44);
end
