function out = smooth_mean2d(A, Nr, Nc)
% SMOOTH_MEAN2D  moving-average smoothing of a 2-D array (used only for the optional field plots).
%   out(i,j) = mean of A over the rectangle of (2*Nr+1) x (2*Nc+1) elements centred on (i,j);
%   at the edges the part of the rectangle that fits is used (edge normalisation).
if nargin < 3, Nc = Nr; end
K = ones(2*Nr+1, 2*Nc+1);
out = conv2(A, K, 'same') ./ conv2(ones(size(A)), K, 'same');
end
