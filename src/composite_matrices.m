function hle = composite_matrices(h)
% COMPOSITE_MATRICES  mass, stiffness and damping matrices of the non-SMA part (GFRP skins)
%   of the GFRP/SMA/GFRP sandwich beam (thesis chapter 4).
%
%   h fields (defaults of the thesis in cases/case_composite.m):
%     N, L, b, rho (GFRP), E (GFRP), thick (total sandwich thickness),
%     thick2 (SMA core thickness), rho2, E2 (SMA, only used to get the Rayleigh coefficients),
%     zeta  (Rayleigh damping ratio imposed on modes 1 and 2).
%   The skins are "total beam minus the SMA-thickness slab"; the SMA slab itself is
%   simulated by the nonlinear solver, the Rayleigh coefficients are computed on the
%   whole sandwich (GFRP skins + SMA core).
N = h.N;
b = h.b;
Lk = h.L/(N-1);

[Mk1,Kk1] = element_matrices(h.rho,  b*h.thick,  h.E,  1/12*b*h.thick^3,  Lk);   % GFRP, full thickness
[Mk2,Kk2] = element_matrices(h.rho,  b*h.thick2, h.E,  1/12*b*h.thick2^3, Lk);   % GFRP, thickness of the SMA core
[Mk3,Kk3] = element_matrices(h.rho2, b*h.thick2, h.E2, 1/12*b*h.thick2^3, Lk);   % SMA core (linear, for Rayleigh only)
Mk = Mk1 - Mk2 + Mk3;
Kk = Kk1 - Kk2 + Kk3;
Mk_vero = Mk1 - Mk2;
Kk_vero = Kk1 - Kk2;

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
Extr_tr = permute(Extr, [2,1,3]);

M = zeros(2*N,2*N,N-1); K = M; M_vero = M; K_vero = M;
for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
    K(:,:,k) = Extr_tr(:,:,k) * Kk * Extr(:,:,k);
end
M = sum(M,3);
K = sum(K,3);
MFF_HLE = M(1:2*N-2,1:2*N-2);
KFF_HLE = K(1:2*N-2,1:2*N-2);

[eigenvectors, eigenvalues]=eig(MFF_HLE\KFF_HLE);
freq=sqrt(diag(eigenvalues))/2/pi;
[dummy,ind]=sort(freq);
h1=h.zeta; h2=h.zeta;
hh=[h1;h2];
om1=2*pi*dummy(1);
om2=2*pi*dummy(2);
H=[1/(2*om1) om1/2;1/(2*om2) om2/2];
coef=H\hh;
alpha=coef(1);
beta=coef(2);

for k = 1:N-1
    M_vero(:,:,k) = Extr_tr(:,:,k) * Mk_vero * Extr(:,:,k);
    K_vero(:,:,k) = Extr_tr(:,:,k) * Kk_vero * Extr(:,:,k);
end
M_vero = sum(M_vero,3);
K_vero = sum(K_vero,3);

MFF_HLE = M_vero(1:2*N-2,1:2*N-2);
KFF_HLE = K_vero(1:2*N-2,1:2*N-2);
CFF_HLE=alpha*MFF_HLE+beta*KFF_HLE;
hle = struct('MFF_HLE',MFF_HLE,'KFF_HLE',KFF_HLE,'CFF_HLE',CFF_HLE,'alpha',alpha,'beta',beta);
end

function [Mk,Kk] = element_matrices(rho,Area,E,J,Lk)
Mk = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
Kk = E*J * [  +12./Lk^3  6./Lk^2     -12./Lk^3  6./Lk^2
               6./Lk^2  +4./Lk       -6./Lk^2  +2./Lk
               -12./Lk^3  -6./Lk^2   +12./Lk^3  -6./Lk^2
               6./Lk^2  +2./Lk       -6./Lk^2  +4./Lk ] ;
end
