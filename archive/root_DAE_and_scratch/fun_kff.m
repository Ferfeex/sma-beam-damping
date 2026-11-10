function KFF = fun_kff(sigma,eps,xisp,xism,xiT,Extr_tr,Extr,N,x,y,Ea,Emp,Emm,EmT,eps_Lp,eps_Lm,b,ssp,sfp,ssm,sfm,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0)
global sigma_old
% DXI ARE THE DERIVATIVES OF MARTENSITIC FRACTIONS XI, FOR TANGENT MATRIX
indp = find(sigma>ssp & sigma < sfp & sigma >= sigma_old);
indm = find(sigma<ssm & sigma > sfm & sigma <= sigma_old);
Dxisp = zeros(750,1); Dxism = zeros(750,1); DxiT = zeros(750,1);
if any(indp)~=0
    Dxisp(indp) = pi / (ssp-sfp) * sin(pi*(sigma(indp)-ssp)/(sfp-ssp)) .* 1/2.*(xisp0(indp)-1);
    Dxism(indp) = xism0(indp) ./ (xisp0(indp) - 1) .* Dxisp(indp);
    DxiT(indp) = xiT0(indp) ./ (xisp0(indp) - 1) .* Dxisp(indp);
end
if any(indm)~=0
    Dxism(indm) = pi / (ssm-sfm) * sin(pi*(sigma(indm)-ssm)/(sfm-ssm)) .* 1/2.*(xism0(indm)-1);
    Dxisp(indm) = xisp0(indm) ./ (xism0(indm) - 1) .* Dxism(indm);
    DxiT(indm) = xiT0(indm) ./ (xism0(indm) - 1) .* Dxism(indm);
end
% YOUNG'S AND ITS DERIVATIVE
E = Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea);
dE = Dxisp*Emp + Dxism*Emm + DxiT*EmT;
% TANGENT MATRIX C
C_tot = E./(1-(eps-eps_Lp*xisp-eps_Lm*xism).*dE+E.*(eps_Lp*Dxisp+eps_Lm*Dxism));
Ck = reshape(C_tot,length(y),length(x),N-1);
% STIFFNESS MATRIX FOR K-TH ELEMENT
Kk = zeros(4,4,N-1);
for k = 1:N-1
    Kk(:,:,k) = b*[trapz(x,trapz(-y,(Ck(:,:,k).*BB11))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB12))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB13))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB14)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB21))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB22))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB23))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB24)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB31))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB32))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB33))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB34)))
        trapz(x,trapz(-y,(Ck(:,:,k).*BB41))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB42))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB43))) , trapz(x,trapz(-y,(Ck(:,:,k).*BB44)))];
end
% ASSEMBLING THE STIFFNESS MATRICES
K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    K(:,:,k) = Extr_tr(:,:,k) * Kk(:,:,k) * Extr(:,:,k);
end
K = sum(K,3);
KFF = K(1:2*N-2,1:2*N-2);
end
