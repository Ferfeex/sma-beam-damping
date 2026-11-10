clear all
clc
close all
% GEOMETRY & FEM
N = 16; % Numero nodi
L = 1000; % [mm] lunghezza trave
b = 30; % b e h cross section
thick = 20;
Area = b * thick;
J = 1/12 * b * thick^3;
rho = 6.1e-09; % ton/mm^3 density
Lk = L/(N-1);
csi = linspace(0,Lk,5)'; % x discretization within single element
y = linspace(10,-10,10)'; % y discretization
% SMA PROPERTIES
Ea = 34.9e3; %[MPa]
Emp = 26e3;
Emm = 33e3;
EmT = 26e3;
ssp = 0; %[MPa]
sfp = 188;
ssm = 0;
sfm = -220;
eps_Lp = 0.016;
eps_Lm = -0.0146;
% ELEMENT MASS MATRIX
Mk = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
% VERTICAL LOAD
F5 = zeros(4,1);
F5(1) = -5;
%ASSEMBLY
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
M = zeros(2*N,2*N,N-1);
for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
end
M = sum(M,3);
F = Extr_tr(:,:,ceil(N/2)) * F5;
MFF = M(1:2*N-2,1:2*N-2);
FF = F(1:2*N-2);
% CALCULATION OF B AND B' CONSTANT MATRICES FOR K NONLINEAR MATRIX
%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS FOR SIMULATION
X0 = zeros(3810,1);
X0(1:30) = 0.01*[196.286220936686,2.90106152335298,383.983151511633,2.70922553770192,554.878722267273,2.39887310003360,701.498509654886,1.98362709781420,817.433353149620,1.48170374163203,897.619677242082,0.915104715279853,938.560225853662,0.308644277012870,938.475665006134,-0.311143813136586,897.380646506695,-0.917171229075329,817.081248547045,-1.48297970406239,701.094274038331,-1.98389224637425,554.492269619767,-2.39808244477330,383.681262795266,-2.70751656912837,196.121022338791,-2.89872987584355,2.96595696909290,-2.96340518254054]';
X0(1561:2310) = 1;
% X0(1:30) NODES VELOCITY AS PERTURBATION (NO FORCE WILL BE INCLUDED)
% X0(1561:2310) XIT --> XISP AND XISM = 0 ; XIT = 1
global sigma_old ex_indp ex_indm xisp_old xism_old xiT_old
sigma_old = X0(2311:3060);  xisp_old = X0(61:810); xism_old = X0(811:1560); xiT_old = X0(1561:2310);
xisp0 = X0(61:810); xism0 = X0(811:1560); xiT0 = X0(1561:2310);
% SIMULATION
T0 = 0; Tf = 0.5;
Mass = zeros(3810);
Mass(1:60,1:60) = eye(60);
load('Sparsity.mat')
options = odeset('Mass',Mass,'MassSingular','yes','OutputFcn',@myfun,'JPattern',Sparsity,'OutputSel',61:3060,'MStateDependence','none','RelTol',1,'AbsTol',0.1);
%%%%%%%%% PRIMA CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tf],X0,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT = T(1:end-1);
XX = X(1:end-1,:);
X0 = XX(end,:)';
T0 = TT(end);
% Tff = T(end);
% options2 = odeset('Mass',Mass,'MassSingular','yes','OutputFcn',@myfun,'JPattern',Sparsity,'OutputSel',61:3060,'MStateDependence','none','RelTol',0.1,'AbsTol',0.01);
%SECONDA CALL SENZA AGGIORNAMENTO XI PER REFINEMENT ULTIMO STEP DELLA PRIMA
% [T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tff],X0,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT = [TT ; T(1:end-1)];
% XX = [XX ; X(1:end-1,:)];
% X0 = XX(end,:)';
% T0 = TT(end);
% Tff = T(end);
% sigma_old = X0(2311:3060);
% options3 = odeset('Mass',Mass,'MassSingular','yes','OutputFcn',@myfun2,'JPattern',Sparsity,'OutputSel',61:3060,'MStateDependence','none','RelTol',0.1,'AbsTol',0.01);
%TERZA CALL SENZA AGGIORNAMENTO XI PER REFINEMENT ULTIMO STEP DELLA PRIMA
% [T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tff],X0,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT = [TT ; T(1:end-1)];
% XX = [XX ; X(1:end-1,:)];
% X0 = XX(end,:)';
% T0 = TT(end);
xisp0(ex_indp) = X0(60+ex_indp); xism0(ex_indp) = X0(810+ex_indp); xiT0(ex_indp) = X0(1560+ex_indp);
xisp0(ex_indm) = X0(60+ex_indm); xism0(ex_indm) = X0(810+ex_indm); xiT0(ex_indm) = X0(1560+ex_indm);
sigma_old = X0(2311:3060);
while TT(end) < Tf
   [T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tf],X0,options);
   if T(end)<Tf
       TT = [TT ; T(1:end-1)];
       XX = [XX ; X(1:end-1,:)];
       X0 = XX(end,:)';
       T0 = TT(end);
%        Tff = T(end);
%        sigma_old = X0(2311:3060);
%        options2 = odeset('Mass',Mass,'MassSingular','yes','OutputFcn',@myfun,'JPattern',Sparsity,'OutputSel',61:3060,'MStateDependence','none','RelTol',0.1,'AbsTol',0.01);
%        [T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tff],X0,options2);
%        TT = [TT;T(1:end-1)];
%        XX = [XX;X(1:end-1,:)];
%        X0 = XX(end,:)';
%        T0 = TT(end);
%        Tff = T(end);
%        sigma_old = X0(2311:3060);
%        options3 = odeset('Mass',Mass,'MassSingular','yes','OutputFcn',@myfun2,'JPattern',Sparsity,'OutputSel',61:3060,'MStateDependence','none','RelTol',0.1,'AbsTol',0.01);
%       if Tff-T0 >=0.001
%           [T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0),[T0 Tff],X0,options);
%           TT = [TT ; T(1:end-1)];
%           XX = [XX ; X(1:end-1,:)];
%           X0 = XX(end,:)';
%           T0 = TT(end);
%       end
       xisp0(ex_indp) = X0(60+ex_indp); xism0(ex_indp) = X0(810+ex_indp); xiT0(ex_indp) = X0(1560+ex_indp);
       xisp0(ex_indm) = X0(60+ex_indm); xism0(ex_indm) = X0(810+ex_indm); xiT0(ex_indm) = X0(1560+ex_indm);
       sigma_old = X0(2311:3060);
   end
end
dy = XX(:,1:30);
y = XX(:,31:60);
xisp = XX(:,61:810);
xism = XX(:,811:1560);
xiT = XX(:,1561:2310);
sigma = XX(:,2311:3060);
eps = XX(:,3061:3810);
figure
plot(TT(1:end-1),y(:,15));
grid on
xlabel('t[s]')
ylabel('displacement middle [mm]')
% STRAIN FIELD & MARTENSITIC FRACTIONS
% meshare(eps,xisp,xism,xiT)
% DISPLACEMENT FIELD
% csi = linspace(0,Lk,30)';
% NW1=2*(csi/Lk).^3-3*(csi/Lk).^2+1;
% NW2=Lk*((csi/Lk).^3-2*(csi/Lk).^2+csi/Lk);
% NW3=-2*(csi/Lk).^3+3*(csi/Lk).^2;
% NW4=Lk*((csi/Lk).^3-(csi/Lk).^2);
% NN=[NW1,NW2,NW3,NW4];
% axis([0,1000,-0.0005,0.0005]);
% h=gca;
% hold on
% h.YLimMode='manual';
% for ind = 1:1:length(X(:,1))
%     XA = X(ind,:)';
%     XF = XA(2*N-1:60);
%     XA = [XF ; zeros(2,1)];
%     for k  = 1:N-1
%         u(:,k) = NN * Extr(:,:,k) * XA;
%     end
%     ww = u(:,:);
%     ww = ww(:);
%     x = linspace(0,L,length(ww))';
%     plot(x,ww,'b',x,zeros(length(x)),'r')
%     grid on
%     legend('deformed','undeformed')
%     title('displacement field')
%     xlabel('length [mm]')
%     ylabel('vertical displacement [mm]')
%     drawnow
%     cla
% end
% hold off
