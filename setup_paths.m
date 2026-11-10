function setup_paths()
% SETUP_PATHS  add the repository folders to the Octave/MATLAB path (called by run_case).
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root,'src'));
addpath(fullfile(root,'cases'));
addpath(fullfile(root,'third_party'));
addpath(fullfile(root,'brinson'));
% GNU Octave has no ode23t: add the stand-in ONLY under Octave so it never shadows MATLAB's.
if exist('OCTAVE_VERSION','builtin') ~= 0
    addpath(fullfile(root,'compat','octave'));
end
end
