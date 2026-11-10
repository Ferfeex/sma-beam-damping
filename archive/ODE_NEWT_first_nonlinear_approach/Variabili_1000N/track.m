clear
clc
% load('X100.mat')
load('matlab.mat')
e1 = zeros(750,1);
xp1 = zeros(750,1);
xm1 = zeros(750,1);
xt1 = ones(750,1);
s1 = zeros(750,1);
xisp0 = zeros(750,1);
xism0 = zeros(750,1);
xiT0 = ones(750,1);
for k = 1:length(X(:,1))
    [epsilon,sigma,xisp,xism,xiT] = evolution(X(k,31:60),xisp0,xism0,xiT0,xp1(:,end),xm1(:,end),xt1(:,end),s1(:,end));
    e1 = [e1 , epsilon];
    xp1 = [xp1 , xisp];
    xm1 = [xm1 , xism];
    xt1 = [xt1 , xiT];
    s1 = [s1 , sigma];
    ds = [zeros(750,1) diff(s1')'];
    sigma_old = s1(:,end);
    prodd = ds(:,end-1).*ds(:,end);
    ind = find(prodd<0);
    if any(ind)
        xisp0(ind,1) = xp1(ind,end-1);
        xism0(ind,1) = xm1(ind,end-1);
        xiT0(ind,1) = xt1(ind,end-1);
        k = k-1;
        e1(:,end) = e1(:,end-1);
        xp1(:,end) = xp1(:,end-1);
        xm1(:,end) = xm1(:,end-1);
        xt1(:,end) = xt1(:,end-1);
        s1(:,end) = s1(:,end-1);
        [epsilon,sigma,xisp,xism,xiT] = evolution(X(k,31:60),xisp0,xism0,xiT0,xp1(:,end),xm1(:,end),xt1(:,end),s1(:,end));
        e1 = [e1 , epsilon];
        xp1 = [xp1 , xisp];
        xm1 = [xm1 , xism];
        xt1 = [xt1 , xiT];
        s1 = [s1 , sigma];
        ds = [zeros(750,1) diff(s1')'];
        k = k+1;
    end
    k
end
meshare(e1,xp1,xm1,xt1)
