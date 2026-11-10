function [xvect,niter] = Newton(x0,maxit,toll,y,dy,m)
niter = 0;
xvect = x0;
err_norm = toll + 1;
while niter <= maxit && err_norm >= toll
    x = x0 - m * y(x0) / dy(x0);
    xvect(end+1,1) = x;
    err_norm = norm(x0 - x) / norm(x0);
    x0 = x;
end
niter = length(xvect);
end
