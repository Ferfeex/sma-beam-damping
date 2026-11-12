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
[same,detail] = compare_runs(o,r); check('base (Tf=0.03) ~= original (rel 1e-5)', same, detail);

% 2b. field figures: plot_fields writes ONE martensite figure and ONE strain/deformation figure (no empty window)
od = tempname(); mkdir(od); nfig0 = numel(findall(0,'type','figure'));
rf = run_case('base','Tf',0.03,'figures','save','plot_fields',true,'save_mat',false,'outdir',od);
check('plot_fields: <case>_fields.png and <case>_strain.png written, no figure left open', ...
    exist(fullfile(od,'base_fields.png'),'file')==2 && exist(fullfile(od,'base_strain.png'),'file')==2 && numel(findall(0,'type','figure'))==nfig0, ...
    sprintf('(%d figures open)', numel(findall(0,'type','figure'))));

% 2c. mid-span value: deflection at exactly x = L/2 (Hermite) - consistent with the strain_figure shape, the
%     node value, and exact for synthetic data (src/midspan_deflection.m; dof map of build_beam.m)
Nn = 16; Ln = 1000; Lk = Ln/(Nn-1); xn = (0:Nn-1)'*Lk;
wn = sin(pi*xn/Ln); thn = pi/Ln*cos(pi*xn/Ln); un = zeros(1,2*Nn-2);
for k = 2:Nn-1, un(2*k-3) = wn(k); un(2*k-2) = thn(k); end
un(2*Nn-3) = thn(1); un(2*Nn-2) = thn(Nn);
[wmid,wnode,minfo] = midspan_deflection([zeros(1,2*Nn-2) un], Nn, Ln);
check('midspan_deflection: sin(pi x/L) -> 1 at x = L/2 (Hermite, N=16), node 9 is lower', abs(wmid-1) < 1e-4 && abs(wnode-sin(pi*xn(9)/Ln)) < 1e-12 && minfo.x_node > 500, ...
    sprintf('w(L/2) = %.7f, node %d (x = %.1f mm) = %.5f', wmid, minfo.node, minfo.x_node, wnode));
Nn = 15; Lk = Ln/(Nn-1); xn = (0:Nn-1)'*Lk; wn = sin(pi*xn/Ln); thn = pi/Ln*cos(pi*xn/Ln); un = zeros(1,2*Nn-2);
for k = 2:Nn-1, un(2*k-3) = wn(k); un(2*k-2) = thn(k); end
un(2*Nn-3) = thn(1); un(2*Nn-2) = thn(Nn);
[wmid,~,minfo] = midspan_deflection([zeros(1,2*Nn-2) un], Nn, Ln);
check('midspan_deflection: N-1 even -> the middle node itself', abs(wmid-1) < 1e-12 && minfo.x_node == 500, sprintf('node %d at x = %g mm', minfo.node, minfo.x_node));
fh = figure('visible','off');
[xs_,ws_] = strain_figure(rf.fields.epsyylon, rf.XX, 16, 1000, 10, rf.TT(end), fh); close(fh);
[~,im_] = min(abs(xs_-500));
check('res.midspan(end) == deformed shape of strain_figure at x = L/2', abs(xs_(im_)-500) < 1e-9 && abs(ws_(im_)-rf.midspan(end)) < 1e-9*max(1,abs(rf.midspan(end))), ...
    sprintf('%.6f vs %.6f mm', ws_(im_), rf.midspan(end)));
mm_ = abs(rf.midspan_node) > 0.5*max(abs(rf.midspan_node));
rat_ = rf.midspan(mm_)./rf.midspan_node(mm_);
check('mode-1 midspan is ~0.5-0.7 % larger than node 9', all(rat_ > 1.003 & rat_ < 1.010) && rf.x_midspan_node > 530 && rf.x_midspan_node < 535, ...
    sprintf('ratio %.4f..%.4f, node 9 at x = %.1f mm', min(rat_), max(rat_), rf.x_midspan_node));

% 3. global state is reset between runs: same case twice in one session -> identical
r2 = run_case('base','Tf',Tf,'figures','none','save_mat',false);
check('re-run in the same session identical', isequal(r.XX,r2.XX), sprintf('max |dXX| = %g', max(abs(r.XX(:)-r2.XX(:)))));

% 4. composite
o = run_original('Loss_Composite_FEM1D.m', {'data/p.mat','p.mat'; 'tests/reference/Matrices_composite_original.mat','Matrices.mat'; fullfile(orig,'Loss_Composite_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf, cell(0,2));
% (1) with the 2020 Matrices.mat the result is bit-identical to the original script
cfgc = case_composite(); cfgc.time.Tf = Tf;
cfgc.hle.matrices_file = fullfile(root,'tests','reference','Matrices_composite_original.mat');
r = run_case(cfgc,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('composite (Tf=0.03, original Matrices.mat) ~= original (rel 1e-5)', same, detail);
% (2) with the matrices recomputed by composite_matrices.m: ~1e-12 round-off in the Rayleigh coefficients
%     (the 2020 file was written by MATLAB), which gives ~1e-9 relative differences in the result
r = run_case('composite','Tf',Tf,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('composite (Tf=0.03, recomputed matrices) ~= original', same, detail);

% 5. forced pulse (Loss/Transient): short run, includes the ramp and the 0.02 s release
Tf5 = 0.025;
o = run_original('Loss_Transient_FEM1D.m', {'data/p.mat','p.mat'; fullfile(orig,'Loss_Transient_myfun.m'),'myfun.m'; fullfile(orig,'Loss_Transient_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf5, cell(0,2));
r = run_case('pulse','Tf',Tf5,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('pulse (Tf=0.025) ~= original (rel 1e-5)', same, detail);

% 6. super-elastic pre-strain (ode15s of the original replaced by the ode23t stand-in, as in the port)
Tf6 = 0.02;
o = run_original('Loss_SE_pre-strain_FEM1D.m', {'data/p_SE.mat','p.mat'; fullfile(orig,'Loss_SE_pre-strain_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf6, {'ode15s(@','ode23t(@'});
r = run_case('prestrain_se','Tf',Tf6,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('prestrain_se (Tf=0.02) ~= original (rel 1e-5)', same, detail);

if strcmp(mode,'full')
    % 7. full base case: peaks vs thesis Fig. 4.3 and the stored reference run of the Octave port
    r = run_case('base','figures','save','save_mat',false);
    pk = r.peaks(:,2)';
    expect = [7.13 5.52 4.17 3.09 2.25 1.62];   % deflection at x = L/2 (Hermite), corrected (||) Newton loop; node 9 (before 13 Nov 2026) gave 7.09 5.49 4.15 3.07 2.24 1.61
    check('base full run: ring-down peaks', numel(pk) >= 6 && all(abs(pk(1:6)-expect) < 0.02), ...
        sprintf('peaks [mm]: %s (expected %s; thesis Fig. 4.3: ~5.5 4.15 3.1 2.2 1.6 after the first)', mat2str(pk(1:min(6,end)),3), mat2str(expect)));
    r = run_case('composite','figures','save','save_mat',false);
    check('composite full run finished', r.TT(end) > 0.7, sprintf('t_end = %.3f s, %d steps, %.0f s', r.TT(end), numel(r.TT), r.runtime_s));
end
ok = (nfail == 0);
if ok, fprintf('\nALL TESTS PASSED\n'); else, fprintf('\n%d TEST(S) FAILED\n', nfail); end
end

function [same,detail] = compare_runs(o,r,tol)
if nargin < 3, tol = 1e-5; end   % relative to max|X|; round-off level (the fast path of 2026 sums in a different order, see PERF_NOTES.md)
% o: original (TT,XX); r: refactored result
if numel(o.TT) ~= numel(r.TT)
    same = false; detail = sprintf('different number of steps: original %d, new %d', numel(o.TT), numel(r.TT)); return
end
dX = max(abs(o.XX(:)-r.XX(:))); dT = max(abs(o.TT(:)-r.TT(:)));
scale = max(abs(o.XX(:)));
same = dX <= tol*scale && dT <= tol*max(abs(o.TT(:)));
detail = sprintf('%d steps, max|dXX| = %.2e (rel %.1e), max|dT| = %.1e', numel(r.TT), dX, dX/scale, dT);
end
