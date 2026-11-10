function eps = fun_eps(B,N,X,Extr)
% DISPLACEMENT FIELD
XA = [X;zeros(2,1)];
X_k = zeros(4*(N-1),1);
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
end
% STRAIN FIELD
eps = B * X_k;
end
