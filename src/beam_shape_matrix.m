function [x, Hm] = beam_shape_matrix(N, L, nx)
% BEAM_SHAPE_MATRIX  linear map from the free displacement dofs u (2N-2 values, = res.XX(:,2N-1:end).')
%   to the deflection w(x) on a grid of nx points (nx odd, so that x = L/2 is on the grid):
%       w = Hm * u      (Hm is nx x (2N-2),  x is nx x 1)
%   Inside each element the deflection is the cubic Hermite interpolation of w and theta = dw/dx of its
%   two nodes (the shape functions of the Euler-Bernoulli element of the FE model, as in
%   src/midspan_deflection.m and src/strain_figure.m). DOF MAP (src/build_beam.m): node k = 2..N-1 has
%   w = u(2k-3), theta = u(2k-2); w_1 = w_N = 0 (supports); theta_1 = u(2N-3), theta_N = u(2N-2).
if nargin < 3, nx = 201; end
if mod(nx,2) == 0, nx = nx+1; end
Lk = L/(N-1);
x  = linspace(0, L, nx).';
iw = zeros(N,1);  it = zeros(N,1);                 % index of w_k, theta_k in u (0 = not a dof: w at the supports)
for k = 2:N-1, iw(k) = 2*k-3; it(k) = 2*k-2; end
it(1) = 2*N-3;  it(N) = 2*N-2;
Hm = zeros(nx, 2*N-2);
for i = 1:nx
    e = min(floor(x(i)/Lk + 1e-9) + 1, N-1);
    s = (x(i) - (e-1)*Lk)/Lk;
    h = [1-3*s^2+2*s^3, Lk*(s-2*s^2+s^3), 3*s^2-2*s^3, Lk*(-s^2+s^3)];   % w_e, th_e, w_e+1, th_e+1
    idx = [iw(e) it(e) iw(e+1) it(e+1)];
    for j = 1:4
        if idx(j) > 0, Hm(i,idx(j)) = Hm(i,idx(j)) + h(j); end
    end
end
end
