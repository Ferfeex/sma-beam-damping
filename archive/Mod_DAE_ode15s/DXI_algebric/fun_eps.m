function eps = fun_eps(B,N,X,Extr)
% xe e ye rappresentano le coordinate del sistema di riferimento locale
% dell'elemento finito
X = [X;zeros(2,1)];
X_k = zeros(4*(N-1),1);
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * X;
end
eps = B * X_k;
