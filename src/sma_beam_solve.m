clear all
clc

% GEOMETRY & FEM
N = 16; % Nodes
L = 1000; % Length [mm]
b = 15; % Width [mm]
thick = 10; % Thickness [mm]
Area = b * thick;
J = 1/12 * b * thick^3; % Inertia
Lk = L/(N-1);
csi = linspace(0,Lk,5)'; % x discretization within single element
y = linspace(5,-5,10)'; % y discretization
% SMA PROPERTIES
Ea = 34.9e3;
Emp = 26e3;
Emm = 33e3;
EmT = 26e3;
ssp = 0;
sfp = 188;
ssm = 0;
sfm = -220;
eps_Lp = 0.016;
eps_Lm = -0.0146;
rho = 6.1e-09; % Density [ton/mm^3]
% ELEMENT MASS
Mk = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];
E = EmT;
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
M = zeros(32);
for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
    K(:,:,k) = Extr_tr(:,:,k) * Kk * Extr(:,:,k);
end
M = sum(M,3);
K = sum(K,3);
% FREE DOFS
global KFF CFF LF
MFF = M(1:2*N-2,1:2*N-2);
KFF = K(1:2*N-2,1:2*N-2);
load('p');
LF_max = polyval(p,0.0006);
CFF = KFF * LF_max / (2*pi*10);
LF = LF_max * ones(750,1);
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
% TIME INTERVAL RANGE
T0 = 0; Tf = 0.8; %[s]
% INITIAL CONDITIONS, NODES VELOCITY ~= 0
X0 = 0.46*[196.286220936686,2.90106152335298,383.983151511633,2.70922553770192,554.878722267273,2.39887310003360,701.498509654886,1.98362709781420,817.433353149620,1.48170374163203,897.619677242082,0.915104715279853,938.560225853662,0.308644277012870,938.475665006134,-0.311143813136586,897.380646506695,-0.917171229075329,817.081248547045,-1.48297970406239,701.094274038331,-1.98389224637425,554.492269619767,-2.39808244477330,383.681262795266,-2.70751656912837,196.121022338791,-2.89872987584355,2.96595696909290,-2.96340518254054]';
X0 = [X0;zeros(30,1)];
% INITIAL CONDITIONS VARIABLES
global xisp xism xiT xisp0 xism0 xiT0 indp_old indm_old sigma_old epsyylon xysp xysm xyT s
sigma_old = zeros(750,1); xisp0 = zeros(750,1); s = sigma_old;
xism0 = zeros(750,1); xiT0 = ones(750,1); epsilon = zeros(750,1);
xisp = xisp0; xism = xism0; xiT = xiT0; indp_old = []; indm_old = [];
xysp = xisp; xysm = xism; xyT = xiT; epsyylon = zeros(750,1);
% OPTION FOR ODE SOLVER
options = odeset('OutputFcn',@(T,X,flag)myfun(T,X,flag,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,Extr,Extr_tr,csi,y,B,N,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,p));
% TIME INTEGRATION
[T,X] = ode23t(@(T,X)fun_ode_1D(T,X,MFF,N),[T0 Tf],X0,options);
TT = T;
XX = X;
T0 = T(end);
X0 = X(end,:)';
while T(end) < Tf-0.0001
    [T,X] = ode23t(@(T,X)fun_ode_1D(T,X,MFF,N),[T0 Tf],X0,options);
    TT = [TT ; T(2:end)];
    XX = [XX ; X(2:end,:)];
    T0 = TT(end);
    X0 = X(end,:)';
end
% MID-SPAN VERTICAL DISPLACEMENT
indice = find(XX(:,2*N-2+2*N/2-1) == max(XX(:,2*N-2+2*N/2-1)));
TT = TT(indice:end)-TT(indice);
XX = XX(indice:end,:);
epsyylon = epsyylon(:,indice:end);
xysp = xysp(:,indice:end);
xysm = xysm(:,indice:end);
xyT = xyT(:,indice:end);
s = s(:,indice:end);
plot(TT,XX(:,2*N-2+2*N/2-1),'Linewidth',1.5)
grid on
xlabel('Time t[s]')
ylabel('Displacement \delta [mm]')
title('Mid-span Vertical Displacement')
% meshare(epsyylon,xysp,xysm,xyT,TT)
