function X0 = initial_velocity(N, scale)
% INITIAL_VELOCITY  initial state X0 = [velocities; displacements] of the free vibration.
%   Velocity field = half sine (first bending mode shape) for w and its slope for theta;
%   displacements zero. "scale" multiplies the field (0.46 base, 0.663 composite, 1.485 ...).
%   N = 16 (the thesis case) uses the literal vector of the thesis case.
%   Other N use the same half-sine formula (A = 944.07 mm/s), UNTESTED.
if N == 16
    V = [196.286220936686,2.90106152335298,383.983151511633,2.70922553770192,554.878722267273,2.39887310003360,701.498509654886,1.98362709781420,817.433353149620,1.48170374163203,897.619677242082,0.915104715279853,938.560225853662,0.308644277012870,938.475665006134,-0.311143813136586,897.380646506695,-0.917171229075329,817.081248547045,-1.48297970406239,701.094274038331,-1.98389224637425,554.492269619767,-2.39808244477330,383.681262795266,-2.70751656912837,196.121022338791,-2.89872987584355,2.96595696909290,-2.96340518254054]';
else
    L = 1000; A = 944.07; xk = (1:N-2)*L/(N-1);          % free nodes 2..N-1
    V = zeros(2*N-2,1);
    V(1:2:2*N-4) = A*sin(pi*xk/L);                        % w
    V(2:2:2*N-4) = A*(pi/L)*cos(pi*xk/L);                 % theta (internal nodes)
    V(2*N-3) = A*pi/L;  V(2*N-2) = -A*pi/L;               % theta at the two supports
end
X0 = scale*V;
X0 = [X0; zeros(2*N-2,1)];
end
