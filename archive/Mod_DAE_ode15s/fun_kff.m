function KFF = fun_kff(eps,xi_sp,xi_sm,xi_T,Extr_tr,Extr,N,x,y,Lk,Ea,Emp,Emm,EmT,eps_Lp,eps_Lm,b)
global Dxisp Dxism DxiT
[xe,ye] = meshgrid(x,y); % Il vettore x in input è colonna

BB11 = ye.^2 .* (6*Lk - 12*xe).^2 ./ Lk^6;
BB12 = ye.^2 .* (6*Lk - 12*xe) .*(4*Lk - 6*xe) ./Lk^5;
BB13 = ye.^2 .* (6*Lk - 12*xe) .*(12*xe - 6*Lk) ./Lk^6;
BB14 = ye.^2 .* (6*Lk - 12*xe) .*(2*Lk - 6*xe) ./Lk^5;
BB21 = BB12;
BB22 = ye.^2 .* (4*Lk - 6*xe).^2 ./ Lk^4;
BB23 = -BB12;
BB24 = ye.^2 .* (4*Lk - 6*xe) .* (2*Lk - 6*xe) ./Lk^4;
BB31 = BB13;
BB32 = BB23;
BB33 = BB11;
BB34 = -BB14;
BB41 = BB14;
BB42 = BB24;
BB43 = BB34;
BB44 = ye.^2 .* (2*Lk - 6*xe).^2 ./ Lk^4;

E = Ea + xi_sp * (Emp - Ea) + xi_sm * (Emm - Ea) + xi_T * (EmT - Ea);
dE = Dxisp * (Emp - Ea) + Dxism * (Emm - Ea) + DxiT * (EmT - Ea);

C_tot = E./(1-(eps-eps_Lp*xi_sp-eps_Lm*xi_sm).*dE+E.*(eps_Lp*Dxisp+eps_Lm*Dxism));
Ck = reshape(C_tot,length(y),length(x),N-1);
Kk = zeros(4,4,N-1);
for k = 1:N-1
    Kk(:,:,k) = b*[trapz(x,trapz(-y,(Ck(:,:,k).*BB11))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB12))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB13))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB14)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB21))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB22))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB23))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB24)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB31))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB32))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB33))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB34)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB41))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB42))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB43))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB44)))];
end

K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    K(:,:,k) = Extr_tr(:,:,k) * Kk(:,:,k) * Extr(:,:,k);
end
K = sum(K,3);
KFF = K(1:2*N-2,1:2*N-2);
end
