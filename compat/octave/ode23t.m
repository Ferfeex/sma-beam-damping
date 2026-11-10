function varargout = ode23t(varargin)
% ODE23T  stand-in for MATLAB's ode23t (trapezoidal rule, stiff) in GNU Octave,
% which has no ode23t. Forwards every call to ode23s (Rosenbrock, stiff, supports
% OutputFcn/status restarts). The solver can be changed with the global
% ODE23T_SUBSTITUTE (e.g. 'ode45'). Results are visually the same as in MATLAB
% but not bit-identical (different integrator).
% This folder is added to the path ONLY under Octave (see setup_paths.m), so it
% never shadows MATLAB's own ode23t.
  global ODE23T_SUBSTITUTE
  if isempty(ODE23T_SUBSTITUTE), ODE23T_SUBSTITUTE = 'ode23s'; end
  [varargout{1:max(nargout,1)}] = feval(ODE23T_SUBSTITUTE, varargin{:});
end
