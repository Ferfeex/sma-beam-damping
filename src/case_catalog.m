function C = case_catalog()
% CASE_CATALOG  list of the cases offered by the start window (sma_gui.m) with a description
% and a runtime estimate. Returns a struct array with the fields
%   name       case name as used by run_case('<name>')
%   label      one-line title
%   text       cell array of description lines
%   t_ref      estimated runtime of the full default duration [s], GNU Octave, 1 core
%   Tf_def     default duration of the case [s]
%   slow       true if a default run takes more than ~20 min (none at the moment since the 2026 speed-up)
%   model      'two_phase' | 'one_phase_SE' | 'linear'
%   p_file     file in data/ with the loss-factor polynomial ('' for linear cases)
%   x_ref      damping.x_ref of the case (strain, or stress for the SE case)
% The runtimes are the ones measured in Octave on a Linux box with the fast code path (README, "The cases", PERF_NOTES.md); a laptop will differ.
rows = {
 'base', 'Base case: free vibration of the SMA beam', 10, false, {
   'Cu-Al-Zn-Mn beam 1000 x 15 x 10 mm, simply supported, ferro-elastic.'
   'Released with an initial half-sine velocity; it rings down because every'
   'stress reversal dissipates energy (loss factor eta ~ 0.14, thesis Fig. 4.3).'
   'Expected mid-span peaks: 7.09, 5.49, 4.15, 3.07, 2.24, 1.61 mm.'}
 'composite', 'Composite: GFRP / SMA / GFRP sandwich beam', 2, false, {
   'SMA core 15 x 10 mm between two linear GFRP skins (thesis chapter 4.1,'
   'configuration 1). Free vibration; rings down more slowly than the pure SMA beam.'
   'Expected mid-span peaks: 7.05, 6.25, 5.50, 4.81, 4.19, 3.63 mm.'}
 'hle_linear', 'Linear reference beam (E = 26 GPa, Rayleigh 1.3 %)', 2, false, {
   'The same beam as the base case but linear elastic with Rayleigh damping.'
   'Shows what the SMA damping adds. No loss-factor polynomial is used.'}
 'pulse', 'Force pulse on the SMA beam (200 N / 0.02 s)', 115, false, {
   'Rectangular force pulse on the 30 x 20 mm SMA beam, starting at rest'
   '(thesis 3.3.2, Fig. 3.23). About 2 minutes in Octave with the fast code path'
   '(it took hours before the 2026 speed-up; 64000 solver steps).'}
 'prestrain_se', 'Pre-stressed super-elastic beam (60 MPa)', 45, false, {
   'SMA beam 30 x 20 mm above Af, pre-stressed with 60 MPa, one martensite'
   'fraction (thesis 3.3.3). The loss factor is a function of STRESS (p_SE.mat).'
   'Under a minute in Octave with the fast code path (hours before the 2026 speed-up).'}
 'linear_pulse_incremental', 'Linear GFRP beam, ramp load pulse', 17, false, {
   'Linear GFRP beam (E = 46 GPa), force ramps up to 50 N at t = 0.033 s and is'
   'then released (thesis 3.2, Fig. 3.21). The force term is inferred (KNOWN_ISSUES).'}
 'linear_pulse_rectangular', 'Linear GFRP beam, rectangular load pulse', 17, false, {
   'Linear GFRP beam, rectangular pulse of 37.5 N until t = 0.022 s'
   '(thesis 3.2, Fig. 3.23). The force term is inferred (KNOWN_ISSUES).'}
 };
C = struct('name',{},'label',{},'text',{},'t_ref',{},'Tf_def',{},'slow',{},'model',{},'p_file',{},'x_ref',{});
for k = 1:size(rows,1)
    cfg = feval(['case_' rows{k,1}]);
    C(k).name = rows{k,1};
    C(k).label = rows{k,2};
    C(k).text = rows{k,5};
    C(k).t_ref = rows{k,3};
    C(k).slow = rows{k,4};
    C(k).Tf_def = cfg.time.Tf;
    C(k).model = cfg.model;
    if strcmp(cfg.model,'linear'), C(k).p_file = ''; else, C(k).p_file = cfg.p_file; end
    C(k).x_ref = cfg.damping.x_ref;
end
end
