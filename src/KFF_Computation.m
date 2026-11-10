% Calculation initial KFF
clear all
clc
N = 16;
L = 1000;
b = 30;
Lk = L/(N-1);
csi = linspace(0,Lk,5)';
y = linspace(10,-10,10)';
E = 26e3 * ones(750,1);
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
[CSI,Y] = meshgrid(csi,y);
BB11 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).^2; %Y.^2 .* (6*Lk - 12*CSI).^2 ./ Lk^6;
BB12 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).*(-Y .* (6*CSI/Lk^2 - 4/Lk)); %Y.^2 .* (6*Lk - 12*CSI) .*(4*Lk - 6*CSI) ./Lk^5;
BB13 = -BB11; %Y.^2 .* (6*Lk - 12*CSI) .*(12*CSI - 6*Lk) ./Lk^6;
BB14 = (-Y .* (12*CSI/Lk^3-6/Lk^2)).*(-Y .* (6*CSI/Lk^2 - 2/Lk)); %Y.^2 .* (6*Lk - 12*CSI) .*(2*Lk - 6*CSI) ./Lk^5;
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
C_tot = E;
Ck = reshape(C_tot,length(y),length(csi),N-1);
Kk = zeros(4,4,N-1);
for k = 1:N-1
    Kk(:,:,k) = b*[trapz(csi,trapz(-y,(Ck(:,:,k).*BB11))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB12))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB13))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB14)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB21))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB22))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB23))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB24)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB31))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB32))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB33))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB34)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB41))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB42))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB43))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB44)))];
end
% STIFFNESS MATRIX ASSEMBLY
K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    K(:,:,k) = Extr_tr(:,:,k) * Kk(:,:,k) * Extr(:,:,k);
end
K = sum(K,3);
KFF = K(1:2*N-2,1:2*N-2);
