clear all
clc
close all

N = 16;
rho = 1.88e-9;
E = 46e3;
b = 15;
thick = 12;
Area = b * thick;
J = 1/12 * b * thick^3;
Lk = 1000/15;

Mk1 = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
                   
 
Kk1 = E*J * [  +12./Lk^3  6./Lk^2     -12./Lk^3  6./Lk^2
               6./Lk^2  +4./Lk       -6./Lk^2  +2./Lk 
               -12./Lk^3  -6./Lk^2   +12./Lk^3  -6./Lk^2
               6./Lk^2  +2./Lk       -6./Lk^2  +4./Lk ] ;

           
thick2 = 10;
Area2 = b * thick2;
J2 = 1/12 * b * thick2^3;
           
Mk2 = rho * Area2 * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
                   
 
Kk2 = E*J2 * [  +12./Lk^3  6./Lk^2     -12./Lk^3  6./Lk^2
               6./Lk^2  +4./Lk       -6./Lk^2  +2./Lk 
               -12./Lk^3  -6./Lk^2   +12./Lk^3  -6./Lk^2
               6./Lk^2  +2./Lk       -6./Lk^2  +4./Lk ] ;
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho2 = 6.1e-09;
E2 = 26e3;
Mk3 = rho2 * Area2 * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
                   
 
Kk3 = E2*J2 * [  +12./Lk^3  6./Lk^2     -12./Lk^3  6./Lk^2
               6./Lk^2  +4./Lk       -6./Lk^2  +2./Lk 
               -12./Lk^3  -6./Lk^2   +12./Lk^3  -6./Lk^2
               6./Lk^2  +2./Lk       -6./Lk^2  +4./Lk ] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
h1=0.013; h2=0.013;
h=[h1;h2];
om1=2*pi*dummy(1);
om2=2*pi*dummy(2);
H=[1/(2*om1) om1/2;1/(2*om2) om2/2];
coef=H\h;
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

