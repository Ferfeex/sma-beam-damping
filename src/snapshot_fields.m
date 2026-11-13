function [M, x, t] = snapshot_fields(res, cfg, name, mode)
% SNAPSHOT_FIELDS  martensite fraction at the sampled times as a function of position along the beam.
%   [M,x,t] = snapshot_fields(res, cfg, name, mode)
%     name : 'xisp' (xi_s+, tension), 'xism' (xi_s-, compression) or 'xiT' (twinned martensite)
%     mode : how the thickness is collapsed: 'mean' (average over the section, default), 'top' / 'bottom'
%            (outer fibre y = +h/2 / -h/2), 'max' (largest value in the section)
%   M : nsnap x nc matrix (rows = times t, columns = positions x along the beam, nc = (N-1)*points per element)
%   Uses res.snapshots, recorded when cfg.post.n_snapshots > 0 (run_case does this for plot_fields / animate_fields).
if nargin < 4 || isempty(mode), mode = 'mean'; end
S = res.snapshots;
ny = S.n_y;  nc = size(S.(name),1)/ny;
t = S.t(:);  x = linspace(0, cfg.geom.L, nc);
M = zeros(numel(t), nc);
for j = 1:numel(t)
    F = reshape(S.(name)(:,j), ny, nc);          % rows: y from +h/2 to -h/2 (as in meshare.m)
    switch mode
        case 'mean',   M(j,:) = mean(F,1);
        case 'top',    M(j,:) = F(1,:);
        case 'bottom', M(j,:) = F(end,:);
        case 'max',    M(j,:) = max(F,[],1);
        otherwise, error('snapshot_fields: unknown mode "%s"', mode);
    end
end
end
