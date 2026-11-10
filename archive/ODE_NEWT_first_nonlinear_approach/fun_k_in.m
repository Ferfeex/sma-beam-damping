function [KFF,CFF] = fun_k_in(X)
global xisp0 xism0 xiT0 sigma_old sigma_i xisp xism xiT epsilon xysp xysm xyT
%Material
Ea = 34.9e3;
Emp = 130e3;
Emm = 100e3;
EmT = 70e3;
ssp = 0;
sfp = 120;
ssm = 0;
sfm = -100;
eps_Lp = 0.0012;
eps_Lm = -0.0014;
b = 30;

% Initial displacement field (and null velocity field)
XA = [X(31:60);zeros(2,1)];
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
epsilon = B * X_k;

sigma_i = EmT .* epsilon;
xisp = xysp(:,end);
xism = xysm(:,end);
xiT = xyT(:,end);

index = 1;
indp = find(sigma_i > ssp & sigma_i < sfp & sigma_i >= sigma_old);
indm = find(sigma_i < ssm & sigma_i > sfm & sigma_i <= sigma_old);
i = 0;
Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
while ~isempty(index) 
%     if i < 300000
        xisp(indp) = 1/2 * (1-cos((sigma_i(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos((sigma_i(indp)-ssp)/(sfp-ssp)));
        xiT(indp) = xiT0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
        xism(indp) = xism0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
        xism(indm) = 1/2 * (1-cos((sigma_i(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos((sigma_i(indm)-ssm)/(sfm-ssm)));
        xiT(indm) = xiT0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
        xisp(indm) = xisp0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
        sigmas = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
        
        Dxisp(indp) = 1/2 * pi/(sfp-ssp) * sin(pi*(sigma_i(indp)-ssp)/(sfp-ssp));
        DxiT(indp) = Dxisp(indp) .* (xiT0(indp)./(xisp0(indp)-1));
        Dxism(indp) = Dxisp(indp) .* (xism0(indp)./(xisp0(indp)-1));
        Dxism(indm) = 1/2 * pi/(sfm-ssm) * sin(pi*(sigma_i(indm)-ssm)/(sfm-ssm));
        DxiT(indm) = Dxism(indm) .* (xiT0(indm)./(xism0(indm)-1));
        Dxisp(indm) = Dxism(indm) .* (xisp0(indm)./(xism0(indm)-1));

        dE = Dxisp*(Emp-Ea) + Dxism*(Emm-Ea) + DxiT*(EmT-Ea);
        
        F = sigma_i-sigmas;
        index = find(abs(F)>10^-6);

        dF = 1 - (dE.*(epsilon-eps_Lp*xisp-eps_Lm*xism) + (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)).*(-eps_Lp*Dxisp-eps_Lm*Dxism));
        sigma_i(index) = sigma_i(index) - F(index)./dF(index);
        
        indp = indp(ismember(indp,index));
        indm = indm(ismember(indm,index));
        
        i = i+1;
%     else %Bisection --> SBAGLIATO, NON AGGIORNA DXI CHE SERVE DOPO PER K
%         sigma2 = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
%         xisp2 = xisp;
%         xism2 = xism;
%         xiT2 = xiT;
%         xisp2(indp) = 1/2 * (1-cos((sigma2(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos((sigma2(indp)-ssp)/(sfp-ssp)));
%         xiT2(indp) = xiT0(indp) .* (1 - xisp2(indp))./(1 - xisp0(indp));
%         xism2(indp) = xism0(indp) .* (1 - xisp2(indp))./(1 - xisp0(indp));
%         xism2(indm) = 1/2 * (1-cos((sigma2(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos((sigma2(indm)-ssm)/(sfm-ssm)));
%         xiT2(indm) = xiT0(indm) .* (1 - xism2(indm))./(1 - xism0(indm));
%         xisp2(indm) = xisp0(indm) .* (1 - xism2(indm))./(1 - xism0(indm));
%         sigmas2 = (Ea + xisp2*(Emp-Ea) + xism2*(Emm-Ea) + xiT2*(EmT-Ea)) .* (epsilon - eps_Lp*xisp2 - eps_Lm*xism2);
%         F2 = sigma2 - sigmas2;
%         
%         xisp(indp) = 1/2 * (1-cos((sigma_i(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos((sigma_i(indp)-ssp)/(sfp-ssp)));
%         xiT(indp) = xiT0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
%         xism(indp) = xism0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
%         xism(indm) = 1/2 * (1-cos((sigma_i(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos((sigma_i(indm)-ssm)/(sfm-ssm)));
%         xiT(indm) = xiT0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
%         xisp(indm) = xisp0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
%         sigmas = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
%         F1 = sigma_i-sigmas;
%         
%         aa = F1 .* F2;
%         bb = min(abs(F1),abs(F2));
%         indice1 = find(bb>10^-6 & aa>0);
%         indice2 = find(bb>10^-6 & aa<0);
%         sigma2(indice2) = sigma_i(indice2) + (sigma2(indice2) - sigma_i(indice2))/2;
%         sigma_i(indice1) = sigma2(indice1) + 2*(sigma2(indice1) - sigma_i(indice1));
%         indtol = find(bb<10^-6 & abs(F1)>abs(F2));
%         sigma_i = sigmas;
%         sigma_i(indtol) = sigmas2(indtol);
%         xisp(indtol) = xisp2(indtol);
%         xism(indtol) = xism2(indtol);
%         xiT(indtol) = xiT2(indtol);
%         index = [indice1 ; indice2];
%         indp = indp(ismember(indp,index));
%         indm = indm(ismember(indm,index));
%         i = i+1;
%     end
end

 % Il vettore x in input è colonna
BB11 = Y.^2 .* (6*Lk - 12*CSI).^2 ./ Lk^6;
BB12 = Y.^2 .* (6*Lk - 12*CSI) .*(4*Lk - 6*CSI) ./Lk^5;
BB13 = Y.^2 .* (6*Lk - 12*CSI) .*(12*CSI - 6*Lk) ./Lk^6;
BB14 = Y.^2 .* (6*Lk - 12*CSI) .*(2*Lk - 6*CSI) ./Lk^5;
BB21 = BB12;
BB22 = Y.^2 .* (4*Lk - 6*CSI).^2 ./ Lk^4;
BB23 = -BB12;
BB24 = Y.^2 .* (4*Lk - 6*CSI) .* (2*Lk - 6*CSI) ./Lk^4;
BB31 = BB13;
BB32 = BB23;
BB33 = BB11;
BB34 = -BB14;
BB41 = BB14;
BB42 = BB24;
BB43 = BB34;
BB44 = Y.^2 .* (2*Lk - 6*CSI).^2 ./ Lk^4;

E = Ea + xisp * (Emp - Ea) + xism * (Emm - Ea) + xiT * (EmT - Ea);

C_tot = E./(1-(epsilon-eps_Lp*xisp-eps_Lm*xism).*dE+E.*(eps_Lp*Dxisp+eps_Lm*Dxism));
Ck = reshape(C_tot,length(y),length(csi),N-1);
Kk = zeros(4,4,N-1);
for k = 1:N-1
    Kk(:,:,k) = b*[trapz(csi,trapz(-y,(Ck(:,:,k).*BB11))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB12))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB13))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB14)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB21))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB22))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB23))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB24)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB31))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB32))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB33))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB34)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB41))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB42))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB43))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB44)))];
end

K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    K(:,:,k) = Extr_tr(:,:,k) * Kk(:,:,k) * Extr(:,:,k);
end
K = sum(K,3);
KFF = K(1:2*N-2,1:2*N-2);
CFF = 0.1/100*KFF;
end


