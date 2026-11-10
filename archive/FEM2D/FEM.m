clear
clc
close all

N = 20; % Numero nodi
L = 1000; % [mm] lunghezza trave
b = 30; % b e h cross section
thick = 20;
Area = b * thick;
J = 1/12 * b * thick^3;
rho = 6.1e-09; % ton/mm^3 density
x = linspace(0, L, 500)';
Lk = L/(N-1);
h = 0 : Lk : L;
E = 70000; %MPa

% Structure plot
figure
plot(x,zeros(length(x)), h, zeros(length(h)), 'o')
grid on
xlabel('x[mm]')
ylabel('y[mm]')

% Element matrices
Mk = rho * Area * Lk * [1./3.  0.          0.          1./6.  0.           0.
                    0.     13./35.     11.*Lk/210. 0.     9./70.       -13*Lk/420.
                    0.     11.*Lk/210. Lk^2/105.   0.     13*Lk/420.   -Lk^2/140.
                    1./6.  0.          0.          1./3.  0.           0.
                    0.     9./70.      13*Lk/420.  0.     13./35.      -11.*Lk/210.
                    0.     -13*Lk/420. -Lk^2/140.  0.     -11.*Lk/210. Lk^2/105.   ];

Kk_U = E*Area/Lk * [ 1 0 0 -1 0 0
                0 0 0  0 0 0 
                0 0 0  0 0 0 
                -1 0 0  1 0 0 
                0 0 0  0 0 0 
                0 0 0  0 0 0 ] ;
          
Kk_W = E*J * [ 0.    0.       0.      0.    0.       0.     
               0.  +12./Lk^3  6./Lk^2   0.  -12./Lk^3  6./Lk^2
               0.   6./Lk^2  +4./Lk     0.   -6./Lk^2  +2./Lk
               0.    0.       0.      0.    0.       0. 
               0.  -12./Lk^3  -6./Lk^2  0.  +12./Lk^3  -6./Lk^2
               0.   6./Lk^2  +2./Lk    0.   -6./Lk^2  +4./Lk    ] ;
         
Kk = Kk_U + Kk_W ;

% Vertical load on first node of 10th element
F0 = zeros(6,1);
F5 = zeros(6,1);
F5(2) = -200;

% Assembling the matrices
% Extraction matrices for inner elements
Extr = zeros(6, 3*N, N-1);
for k=2:N-2
    Extr(1:6,3*k-5:3*k,k) = eye(6);
end
% Extraction matrices for boundary elements, hinges on both sides
Extr(4:6,1:3,1) = eye(3);
Extr(3, 3*N - 5, 1) = 1;
Extr(1:2, 3*N-3:3*N-2, 1) = eye(2);

Extr(1:3,3*N-8:3*N-6,N-1) = eye(3);
Extr(6,3*N-4,N-1) = 1;
Extr(4:5,3*N-1:3*N,N-1) = eye(2);

Extr_tr = permute(Extr, [2,1,3]); % Trasposing every extraction matrix

for k = 1:N-1
    M(:,:,k) = Extr_tr(:,:,k) * Mk * Extr(:,:,k);
    K(:,:,k) = Extr_tr(:,:,k) * Kk * Extr(:,:,k);
end

M = sum(M,3);
K = sum(K,3);

F = Extr_tr(:,:,ceil(N/2)) * F5 + Extr_tr(:,:,ceil(N/2)+1) * F5;

MFF = M(1:3*N-4,1:3*N-4);
KFF = K(1:3*N-4,1:3*N-4);
FF = F(1:3*N-4);

% Structural Damping
B=[0.03;0.03];
AA=zeros(2,2);
for ii=1:2
    AA(ii,1)=1/(2* 400 * ii);
    AA(ii,2)= 400 * ii /2;
end
ab=AA\B;
alfa=ab(1);
beta=ab(2);
CFF=alfa*MFF+beta*KFF;
% CFF = zeros(size(CFF));
%
%
% State-space form
A = [-MFF\CFF -MFF\KFF ; eye(3*N-4) zeros(3*N-4)];
B = [inv(MFF) , zeros(3*N-4) ; zeros(3*N-4) , zeros(3*N-4)];
%
T0 = 0; Tf = 1;
X0 = zeros(6*N-8,1);
%
% X0=[0,-0.0782763855484922,-0.00147331552787100,0,-0.154361438431989,-0.00141115010185935,0,-0.226124725337766,-0.00130948071962665,0,-0.291560117149762,-0.00117129571842486,0,-0.348853112106295,-0.00100090017672054,0,-0.396451583075485,-0.000803962120125449,0,-0.433137225281823,-0.000587486508787415,0,-0.458093007412803,-0.000359685508508178,0,-0.470960744072080,-0.000129732674919879,0,-0.471677890620432,0.000104271508488855,0,-0.459883844569467,0.000344038621885447,0,-0.435533024229883,0.000579854110265148,0,-0.399075133736769,0.000802675525394747,0,-0.351410504265104,0.00100446909715868,0,-0.293830515278380,0.00117843653019483,0,-0.227949230970260,0.00131912251037786,0,-0.155632243828101,0.00142241760930341,0,-0.0789274040928020,0.00148549059045145,-0.00149422874025371,0.00150669464673537,0,0.170600136354939,0.00321645660874345,0,0.337261480947940,0.00310420254689339,0,0.496045558083205,0.00291712243859484,0,0.643014520724517,0.00265523115548730,0,0.774231454906348,0.00231854923586744,0,0.885760669131248,0.00190710251910314,0,0.973667959524404,0.00142092158361360,0,1.03402083956896,0.000860040960147391,0,1.06288872240163,0.000224498117885023,0,1.05758647620131,-0.000414792130960655,0,1.02040477558357,-0.000986914595992118,0,0.954879789713548,-0.00149183069043699,0,0.864549664739545,-0.00192950407194126,0,0.752954376459087,-0.00229990171707508,0,0.623635532340633,-0.00260299476628593,0,0.480136135469191,-0.00283875912641010,0,0.326000322822147,-0.00300717584930074,0,0.164773088536230,-0.00310823133272942,0.00325387563193015,-0.00314191740233835];
%
[T,X] = ode15s(@(T,X)fun_ode(T,X,MFF,CFF,KFF,FF,N),[T0 Tf],X0);
plot(T,X(:,3*N-4+3*N/2-1))
grid on
xlabel('t[s]')
ylabel('displacement middle [mm]')

ind = find(X(:,3*N-4+3*N/2-1)==max(X(:,3*N-4+3*N/2-1)));

%% Strain Field
% Displacement for the whole beam
X = X(ind,:)';
XF = X(3*N-3:end);
X = [XF ; zeros(4,1)];
csi = linspace(0,Lk,100)';
NU1=1-(csi)/Lk;
NU2=csi/Lk;
NW1=2*(csi/Lk).^3-3*(csi/Lk).^2+1;
NW2=Lk*((csi/Lk).^3-2*(csi/Lk).^2+csi/Lk);
NW3=-2*(csi/Lk).^3+3*(csi/Lk).^2;
NW4=Lk*((csi/Lk).^3-(csi/Lk).^2);
NN=[NU1,zeros(length(NU1),1),zeros(length(NU1),1),NU2,zeros(length(NU1),1),zeros(length(NU1),1);zeros(length(NU1),1),NW1,NW2,zeros(length(NU1),1),NW3,NW4];
for k  = 1:N-1
    u(:,k) = NN * Extr(:,:,k) * X;
end
uu = u(1:length(csi),:);
uu = uu(:);
ww = u(length(csi)+1:end,:);
ww = ww(:);
x = linspace(0,L,length(uu))';


% % Strain field for the whole beam
% %strain field
y = linspace(10,-10,30)';
B = zeros(length(csi),6,length(y));
B(:,1,:) = -1/Lk * ones(length(csi),length(y));
B(:,4,:) = -B(:,1,:);
dNW1 = 12*csi/Lk^3-6/Lk^2;
B(:,2,:) = dNW1 * (-y');
dNW2 = 6/Lk^2*csi - 4/Lk;
B(:,3,:) = dNW2 * (-y');
dNW3 = -dNW1;
B(:,5,:) = dNW3 * (-y');
dNW4 = 6/Lk^2*csi - 2/Lk;
B(:,6,:) = dNW4 * (-y');
epsilons = zeros(length(y),length(csi),N-1);
for k = 1:N-1
    X_k = Extr(:,:,k) * X;
    for jj = 1:length(y)
        epsilons(jj,:,k) = B(:,:,jj) * X_k;
    end
end
epsilon = epsilons(:,:);

figure
subplot(2,1,1)
plot(x+uu,ww,x,zeros(length(x)))
grid on
legend('deformed','undeformed')
title('displacement field')
xlabel('length [mm]')
ylabel('vertical displacement [mm]')
subplot(2,1,2)
mesh(x,y,epsilon)
xlabel('length [mm]')
ylabel('thickness [mm]')
colorbar
title('Strain field')


%% Non linear part
% First guess of sigma
sigma = E * epsilon;
figure
mesh(x,y,sigma)
xlabel('length [m]')
ylabel('thickness [m]')
colorbar
title('Stress field first guess')
ss = 0;
sf = 50;
Y = cos((sigma - ss)/(sf - ss));
xi = (1-Y)/2;
figure
mesh(x,y,xi)
title('detwinned martensitic fraction field guess')
xlabel('length [mm]')
ylabel('thickness [mm]')
colorbar

index2 = find(xi(1,:) == max(xi(1,:)));
figure
plot(y,xi(:,index2))
title('trend of detwinned martensitic fraction along thickness, beam center')
ylabel('xi')
xlabel('thickness')
