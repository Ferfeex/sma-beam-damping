function ok = dashboard_supported()
% DASHBOARD_SUPPORTED  true if the current graphics toolkit can show the tabbed dashboard (uipanel + uicontrol):
%   MATLAB, or Octave with the Qt toolkit. With gnuplot / FLTK run_case falls back to separate windows.
ok = true;
if exist('OCTAVE_VERSION','builtin') ~= 0
    try, ok = strcmp(graphics_toolkit(),'qt'); catch, ok = false; end
end
end
