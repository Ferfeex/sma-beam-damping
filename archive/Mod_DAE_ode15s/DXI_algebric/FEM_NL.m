clear
clc
close all

N = 16;
L = 1000;
b = 30;
thick = 20;
Area = b * thick;
J = 1/12 * b * thick^3;
rho = 6.1e-09;
Lk = L/(N-1);
Ea = 34.9e3;
Emp = 140e3;
Emm = 90e3;
EmT = 100e3;
ssp = 0;
sfp = 180;
ssm = 0;
sfm = -180;
eps_Lp = 0.0012;
eps_Lm = -0.0014;

% Element matrices
Mk = rho * Area * Lk * [13./35.     11.*Lk/210.    9./70.       -13*Lk/420.
                        11.*Lk/210. Lk^2/105.      13*Lk/420.   -Lk^2/140.
                        9./70.      13*Lk/420.     13./35.      -11.*Lk/210.
                       -13*Lk/420. -Lk^2/140.      -11.*Lk/210. Lk^2/105. ];

% Vertical load on first node of 10th element
F5 = zeros(4,1);
F5(1) = -80;

% Assembling the matrices
% Extraction matrices for inner elements
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
M = zeros(2*N,2*N,N-1);
for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
end
M = sum(M,3);
F = Extr_tr(:,:,ceil(N/2)) * F5;% + Extr_tr(:,:,ceil(N/2)+1) * F5;
MFF = M(1:2*N-2,1:2*N-2);
FF = F(1:2*N-2);

load('CFF.mat')

T0 = 0; Tf = 1;
% X0 = zeros(6060,1);
% X0(1561:2310) = 1;
load('incon.mat')
load('DXI.mat')
X0 = [incon(1:2310);dxi(:);incon(2311:end)];

y = linspace(10,-10,10)';
Mass = zeros(6060);
Mass(1:60,1:60) = eye(60);

csi = linspace(0,Lk,5)';

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

global sigma_old xisp0 xism0 xiT0
sigma_old = X0(2311:3060); xisp0 = X0(61:810); xism0 = X0(811:1560); xiT0 = X0(1561:2310);
 for k = 1:100000
    X0(61:4560) = fun_xi(X0(4561:5310),ssp,sfp,ssm,sfm,X0(2311:3060),X0(3061:3810),X0(3811:4560)); %XI
    X0(4561:5310) = (Ea + X0(61:810) * (Emp - Ea) + X0(811:1560) * (Emm - Ea) + X0(1561:2310) * (EmT - Ea)) .* (X0(5311:6060) - eps_Lp * X0(61:810) - eps_Lm * X0(811:1560)); %SIGMA
 end
load('JPDX.mat')
options = odeset('Events',@reversalEvent,'Mass',Mass,'MassSingular','yes','OutputFcn',@myfun,'JPattern',JPDX);
[T,X] = ode15s(@(T,X)fun_ode_NL(T,X,MFF,CFF,N,FF,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,Lk,b),[T0 Tf],X0,options);

figure
plot(T,X(:,2*N-2+2*N/2-1))
grid on
xlabel('t[s]')
ylabel('displacement middle [mm]')

% ind = find(X(:,2*N-2+2*N/2-1)==max(X(:,2*N-2+2*N/2-1)));
% index = length(X(:,1));
%% Strain Field
% Displacement for the whole beam
% for ind = index
%     XA = X(ind,:)';
%     XF = XA(2*N-1:60);
%     XA = [XF ; zeros(2,1)];
%     csi = linspace(0,Lk,30)';
%     NW1=2*(csi/Lk).^3-3*(csi/Lk).^2+1;
%     NW2=Lk*((csi/Lk).^3-2*(csi/Lk).^2+csi/Lk);
%     NW3=-2*(csi/Lk).^3+3*(csi/Lk).^2;
%     NW4=Lk*((csi/Lk).^3-(csi/Lk).^2);
%     NN=[NW1,NW2,NW3,NW4];
%     for k  = 1:N-1
%         u(:,k) = NN * Extr(:,:,k) * XA;
%     end
%     ww = u(:,:);
%     ww = ww(:);
%     x = linspace(0,L,length(ww))';
%     plot(x,ww,x,zeros(length(x)))
%     grid on
%     legend('deformed','undeformed')
%     title('displacement field')
%     xlabel('length [mm]')
%     ylabel('vertical displacement [mm]')
%     drawnow
% end
% 
% 
% % Strain field for the whole beam
% % %strain field
% y = linspace(10,-10,30)';
% B = zeros(length(csi),4,length(y));
% dNW1 = 12*csi/Lk^3-6/Lk^2;
% B(:,1,:) = dNW1 * (-y');
% dNW2 = 6/Lk^2*csi - 4/Lk;
% B(:,2,:) = dNW2 * (-y');
% dNW3 = -dNW1;
% B(:,3,:) = dNW3 * (-y');
% dNW4 = 6/Lk^2*csi - 2/Lk;
% B(:,4,:) = dNW4 * (-y');
% epsilons = zeros(length(y),length(csi),N-1);
% for k = 1:N-1
%     X_k = Extr(:,:,k) * XA;
%     for jj = 1:length(y)
%         epsilons(jj,:,k) = B(:,:,jj) * X_k;
%     end
% end
% epsilon = epsilons(:,:);
% 
% figure
% mesh(x,y,epsilon)
% xlabel('length [mm]')
% ylabel('thickness [mm]')
% colorbar
% title('Strain field')
