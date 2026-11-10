function [x,niter] = Newton_sempl(x0,maxit,toll,y,dy,m)
niter = 0;
err_norm = toll + 1;
while niter <= maxit && err_norm >= toll
    x = x0 - m * y(x0) / dy(x0);
    err_norm = norm(x0 - x) / norm(x0);
    x0 = x;
    niter = niter +1;
end
end
