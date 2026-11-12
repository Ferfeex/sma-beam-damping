function Mfull = assemble_element_matrix(Cvec, pf)
% ASSEMBLE_ELEMENT_MATRIX  vectorised replacement of the "trapz loops + Extr'*Kk*Extr" blocks of
%   myfun.m / myfun_SE.m.  Cvec: modulus-like quantity at the (N-1)*n_csi*n_y integration points
%   (C_tot for the tangent stiffness, C_tot.*LF/(2*pi*f_ref) for the damping).
%   Element matrix  Kk(a,c) = b * sum_points w(point) * Cvec(point) * BB_ac(point),  the two trapz
%   integrations being the fixed weights w (precomputed in sma_beam_solve -> pf.WB), then the
%   assembly with Extr is a sparse scatter-add (pf.I, pf.J).  Same numbers as the loops up to round-off.
Kall = pf.b * (reshape(Cvec, pf.npt, pf.nel).' * pf.WB);     % nel x 16, column a+4*(c-1) = Kk(a,c)
Kfull = sparse(pf.I, pf.J, reshape(Kall.', [], 1), pf.ndof, pf.ndof);
Mfull = full(Kfull(1:pf.nfree, 1:pf.nfree));
end
