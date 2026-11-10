clear all
clc
close all
% GEOMETRY & FEM
N = 16; % Nodes
L = 1000; % Length [mm]
b = 30; % Width [mm]
thick = 20; % Thickness [mm]
Area = b * thick;
J = 1/12 * b * thick^3; % Inertia
Lk = L/(N-1);
csi = linspace(0,Lk,5)'; % x discretization within single element
y = linspace(10,-10,10)'; % y discretization
% SMA PROPERTIES
Ea = 34.9e3;
Emp = 26e3;
Emm = 33e3;
EmT = 26e3;
ssp = 35;
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
% VERTICAL LOAD, 1ST NODE MID-ELEMENT
F5 = zeros(4,1);
F5(1) = -50;
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
end
M = sum(M,3);
F = Extr_tr(:,:,ceil(N/2)) * F5;
% FREE DOFS
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
% TIME INTERVAL RANGE
T0 = 0; Tf = 0.5; %[s]
% INITIAL CONDITIONS, NODES VELOCITY ~= 0
X0 = [196.286220936686,2.90106152335298,383.983151511633,2.70922553770192,554.878722267273,2.39887310003360,701.498509654886,1.98362709781420,817.433353149620,1.48170374163203,897.619677242082,0.915104715279853,938.560225853662,0.308644277012870,938.475665006134,-0.311143813136586,897.380646506695,-0.917171229075329,817.081248547045,-1.48297970406239,701.094274038331,-1.98389224637425,554.492269619767,-2.39808244477330,383.681262795266,-2.70751656912837,196.121022338791,-2.89872987584355,2.96595696909290,-2.96340518254054]';
X0 = [X0;zeros(30,1)];
% GLOBAL VARIABLES DEFINITION (10 VARIABLES)
global sigma_old differ xisp0 xism0 xiT0 epsilon
global epsyylon xysp xysm xyT ind_old xisp xism xiT s
% GLOBAL VARIABLES INITIALIZATION :
% WITHOUT PRE-STRAIN
sigma_old = zeros(750,1); differ = zeros(750,1); xisp0 = zeros(750,1);
xism0 = zeros(750,1); xiT0 = ones(750,1); epsilon = zeros(750,1);
epsyylon = zeros(750,1); xysp = xisp0; xysm = xism0; xyT = xiT0; ind_old = [];
xisp = xisp0; xism = xism0; xiT = xiT0; s = sigma_old;
% WITH PRESTRAIN
% sigma_old = -8.75*ones(750,1); differ = zeros(750,1); xisp0 = 0.0017*ones(750,1);
% xism0 = zeros(750,1); xiT0 = ones(750,1) - xisp0; epsilon = 1.25e-4*ones(750,1); 
% epsyylon = epsilon; xysp = xisp0; xysm = xism0; xyT = xiT0; ind_old = [];
% xisp = xisp0; xism = xism0; xiT = xiT0; s = sigma_old;
% JACOBIAN SPARSITY PATTERN
JBB = zeros(60);
JBB(1:30,:)=1;
JBB(31:60,1:30)=eye(30);
% OPTION FOR ODE SOLVER
options = odeset('OutputFcn',@myfun,'JPattern',JBB);
% TIME INTEGRATION
[T,X] = ode23t(@(T,X)fun_ode_1D(T,X,MFF,FF,N,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,Extr,Extr_tr,csi,y,B,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44),[T0 Tf],X0,options);
% MID-SPAN VERTICAL DISPLACEMENT
figure
plot(T,X(:,2*N-2+2*N/2-1))
grid on
xlabel('Time t[s]')
ylabel('Displacement \delta [mm]')
title('Mid-span Vertical Displacement')
% BEAM DISPLACEMENT FIELD
% NW1=2*(csi/Lk).^3-3*(csi/Lk).^2+1;
% NW2=Lk*((csi/Lk).^3-2*(csi/Lk).^2+csi/Lk);
% NW3=-2*(csi/Lk).^3+3*(csi/Lk).^2;
% NW4=Lk*((csi/Lk).^3-(csi/Lk).^2);
% NN=[NW1,NW2,NW3,NW4];
% figure
% axis([0,1000,-2,2]);
% h=gca;
% hold on
% h.YLimMode='manual';
% times = 0 : T/20 : T;
% time_steps = find(T==times);
% for ind = time_steps
%     XA = X(ind,:)';
%     XF = XA(2*N-1:end);
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
% % STRAIN FIELD & MARTENSITIC FRACTIONS
meshare(epsyylon,xysp,xysm,xyT)
