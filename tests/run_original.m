function out = run_original(script, support, Tf, replace)
% RUN_ORIGINAL  run an ORIGINAL (2020) script from archive/original_case_scripts in a scratch
% folder, in a child Octave process, with a shortened Tf, and return its TT and XX.
% Used by run_tests.m to prove that the refactored solver gives the same numbers.
%   script   file name in archive/original_case_scripts (copied as FEM1D.m)
%   support  cell array of {file in archive/original_case_scripts, name in scratch}, or
%            {file in data/, name} for p.mat etc. (paths relative to the repo root)
%   Tf       shortened end time
%   replace  cell array of {old,new} text replacements applied to the script (e.g. ode15s -> ode23t)
root = fileparts(fileparts(mfilename('fullpath')));
tmp = tempname(); mkdir(tmp);
txt = fileread(fullfile(root,'archive','original_case_scripts',script));
txt = regexprep(txt,'Tf\s*=\s*0\.8;',sprintf('Tf = %g;',Tf));
for k = 1:size(replace,1), txt = strrep(txt,replace{k,1},replace{k,2}); end
fid = fopen(fullfile(tmp,'FEM1D.m'),'w'); fwrite(fid,txt); fclose(fid);
for k = 1:size(support,1)
    copyfile(fullfile(root,support{k,1}), fullfile(tmp,support{k,2}));
    if strcmp(support{k,2},'myfun.m')
        % The ONLY intended difference between the 2020 original and this code: the Newton loop
        % condition && -> || (KNOWN_ISSUES.md #1). Patch the copy of the original
        % so that the comparison isolates the refactoring (everything else must be bit-identical).
        f = fullfile(tmp,'myfun.m'); t = fileread(f);
        t2 = strrep(t,'while ~isempty(indp) && ~isempty(indm)','while ~isempty(indp) || ~isempty(indm)');
        if strcmp(t,t2), error('run_original: Newton loop line not found in %s', support{k,1}); end
        fid = fopen(f,'w'); fwrite(fid,t2); fclose(fid);
    end
end
fid = fopen(fullfile(tmp,'driver.m'),'w');
fprintf(fid,'set(0,''defaultfigurevisible'',''off'');\naddpath(''%s'');\nFEM1D;\nsave(''-v7'',''orig_out.mat'',''TT'',''XX'');\n', fullfile(root,'compat','octave'));
fclose(fid);
oct = 'octave-cli';
if exist(fullfile(OCTAVE_HOME,'bin','octave-cli'),'file'), oct = fullfile(OCTAVE_HOME,'bin','octave-cli'); end
cmd = sprintf('cd "%s" && "%s" --no-gui --norc driver.m > orig_out.log 2>&1', tmp, oct);
st = system(cmd);
if st ~= 0 || ~exist(fullfile(tmp,'orig_out.mat'),'file')
    error('run_original: child Octave failed, see %s', fullfile(tmp,'orig_out.log'));
end
out = load(fullfile(tmp,'orig_out.mat'));
rmdir(tmp,'s');
end
