function EPS = strain_from_dofs(cfg, U)
% STRAIN_FROM_DOFS  strain at the integration points from nodal displacement dofs, without any field history.
%   EPS = strain_from_dofs(cfg, U)    U: one row per instant, the free displacement dofs (2N-2 columns, i.e.
%   res.XX(:,2N-1:end)).   EPS: (N-1)*n_csi*n_y rows x one column per row of U, the same layout as
%   res.fields.epsyylon (y fastest, then the points along an element, then the elements).
%   epsilon = B*Extr*u (+ sigma0/Ea for the pre-stressed super-elastic case, as in src/myfun_SE.m).
N = cfg.geom.N;
bm = build_beam(cfg, 1);                              % B does not depend on the modulus
Eall = zeros(4*(N-1), 2*N-2);
for k = 1:N-1, Eall(4*(k-1)+1:4*k,:) = bm.Extr(:,1:2*N-2,k); end
EPS = (bm.B * Eall) * U.';
if strcmp(cfg.model,'one_phase_SE'), EPS = EPS + cfg.init.sigma0/cfg.mat.Ea; end
end
