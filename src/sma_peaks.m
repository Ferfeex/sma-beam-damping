function [tp,xp] = sma_peaks(t,x)
% SMA_PEAKS  times and heights of the local maxima of a signal (first sample included
%   if it is a maximum), used to read the ring-down of the mid-span displacement.
x = x(:); t = t(:);
if numel(x) < 3, tp = t(1:0); xp = x(1:0); return, end   % too short to have peaks
loc = find(x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end)) + 1;
if x(1) >= x(2) && x(1) > 0, loc = [1; loc]; end
tp = t(loc); xp = x(loc);
end
