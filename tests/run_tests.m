function ok = run_tests(mode)
% RUN_TESTS  checks that the cleaned-up code still gives the same numbers as the original.
%
%   run_tests            quick tests (~3-6 min): short runs, compared with the ORIGINAL scripts
%   run_tests('full')    + the full 0.8 s base case and composite case (~10 min), peaks vs thesis
%
% What is compared: the refactored solver (run_case / sma_beam_solve) against the unmodified
% original scripts in archive/original_case_scripts/, both run with the same shortened Tf and
% the same Octave ode23t stand-in; the copy of the original is patched ONLY in the Newton loop
% condition (&& -> ||, KNOWN_ISSUES.md #1). Equal to round-off = the numerics of the main method are unchanged.
if nargin < 1, mode = 'quick'; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths(); addpath(fileparts(mfilename('fullpath')));
set(0,'defaultfigurevisible','off');
orig = fullfile('archive','original_case_scripts');
nfail = 0;

    function check(name, cond, detail)
        if cond, fprintf('PASS  %s   %s\n', name, detail); else, fprintf('FAIL  %s   %s\n', name, detail); nfail = nfail + 1; end
    end

% 1. composite matrices vs the ones the author saved by hand in 2020
ref = load(fullfile(root,'tests','reference','Matrices_composite_original.mat'));
cfg = case_composite(); hle = composite_matrices(cfg.hle);
d = max([max(abs(hle.MFF_HLE(:)-ref.MFF_HLE(:))./max(abs(ref.MFF_HLE(:)))), ...
         max(abs(hle.KFF_HLE(:)-ref.KFF_HLE(:))./max(abs(ref.KFF_HLE(:)))), ...
         max(abs(hle.CFF_HLE(:)-ref.CFF_HLE(:))./max(abs(ref.CFF_HLE(:))))]);
check('composite_matrices == original Matrices.mat', d < 1e-9, sprintf('max rel diff %.2e (round-off: the 2020 file was written by MATLAB eig)', d));

% 2. base case: refactored vs original, short run
Tf = 0.03;
o = run_original('Loss_FEM1D_base.m', {'data/p.mat','p.mat'; fullfile(orig,'Loss_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf, cell(0,2));
r = run_case('base','Tf',Tf,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('base (Tf=0.03) == original', same, detail);

% 3. global state is reset between runs: same case twice in one session -> identical
r2 = run_case('base','Tf',Tf,'figures','none','save_mat',false);
check('re-run in the same session identical', isequal(r.XX,r2.XX), sprintf('max |dXX| = %g', max(abs(r.XX(:)-r2.XX(:)))));

% 4. composite
o = run_original('Loss_Composite_FEM1D.m', {'data/p.mat','p.mat'; 'tests/reference/Matrices_composite_original.mat','Matrices.mat'; fullfile(orig,'Loss_Composite_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf, cell(0,2));
% (1) with the 2020 Matrices.mat the result is bit-identical to the original script
cfgc = case_composite(); cfgc.time.Tf = Tf;
cfgc.hle.matrices_file = fullfile(root,'tests','reference','Matrices_composite_original.mat');
r = run_case(cfgc,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('composite (Tf=0.03, original Matrices.mat) == original', same, detail);
% (2) with the matrices recomputed by composite_matrices.m: ~1e-12 round-off in the Rayleigh coefficients
%     (the 2020 file was written by MATLAB), which gives ~1e-9 relative differences in the result
r = run_case('composite','Tf',Tf,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r,1e-7); check('composite (Tf=0.03, recomputed matrices) ~= original', same, detail);

% 5. forced pulse (Loss/Transient): short run, includes the ramp and the 0.02 s release
Tf5 = 0.025;
o = run_original('Loss_Transient_FEM1D.m', {'data/p.mat','p.mat'; fullfile(orig,'Loss_Transient_myfun.m'),'myfun.m'; fullfile(orig,'Loss_Transient_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf5, cell(0,2));
r = run_case('pulse','Tf',Tf5,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('pulse (Tf=0.025) == original', same, detail);

% 6. super-elastic pre-strain (ode15s of the original replaced by the ode23t stand-in, as in the port)
Tf6 = 0.02;
o = run_original('Loss_SE_pre-strain_FEM1D.m', {'data/p_SE.mat','p.mat'; fullfile(orig,'Loss_SE_pre-strain_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf6, {'ode15s(@','ode23t(@'});
r = run_case('prestrain_se','Tf',Tf6,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('prestrain_se (Tf=0.02) == original', same, detail);

if strcmp(mode,'full')
    % 7. full base case: peaks vs thesis Fig. 4.3 and the stored reference run of the Octave port
    r = run_case('base','figures','save','save_mat',false);
    pk = r.peaks(:,2)';
    expect = [7.09 5.49 4.15 3.07 2.24 1.61];   % with the corrected (||) Newton loop; the original && gave 5.494 / 4.155 / 3.076
    check('base full run: ring-down peaks', numel(pk) >= 6 && all(abs(pk(1:6)-expect) < 0.02), ...
        sprintf('peaks [mm]: %s (expected %s; thesis Fig. 4.3: ~5.5 4.15 3.1 2.2 1.6 after the first)', mat2str(pk(1:min(6,end)),3), mat2str(expect)));
    r = run_case('composite','figures','save','save_mat',false);
    check('composite full run finished', r.TT(end) > 0.7, sprintf('t_end = %.3f s, %d steps, %.0f s', r.TT(end), numel(r.TT), r.runtime_s));
end
ok = (nfail == 0);
if ok, fprintf('\nALL TESTS PASSED\n'); else, fprintf('\n%d TEST(S) FAILED\n', nfail); end
end

function [same,detail] = compare_runs(o,r,tol)
if nargin < 3, tol = 1e-9; end
% o: original (TT,XX); r: refactored result
if numel(o.TT) ~= numel(r.TT)
    same = false; detail = sprintf('different number of steps: original %d, new %d', numel(o.TT), numel(r.TT)); return
end
dX = max(abs(o.XX(:)-r.XX(:))); dT = max(abs(o.TT(:)-r.TT(:)));
scale = max(abs(o.XX(:)));
same = dX <= tol*scale && dT <= 1e-9;
detail = sprintf('%d steps, max|dXX| = %.2e (rel %.1e), max|dT| = %.1e', numel(r.TT), dX, dX/scale, dT);
end
