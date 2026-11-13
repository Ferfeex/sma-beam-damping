function s = sma_section_label(cfg, kind)
% SMA_SECTION_LABEL  text that says which part of the section a martensite / strain plot refers to.
%   The martensite fractions and the strain maps are evaluated at the integration points of the SMA section only
%   (cfg.geom.thick, height h = +-thick/2). In the sandwich beam (cfg.hle) the GFRP skins around it are linear
%   elastic: they have no martensite and are not part of these plots.
%     sma_section_label(cfg,'axis')    y-axis label  'h [mm]'  or  'h [mm] (SMA core)'
%     sma_section_label(cfg,'thickness')  the same with the word 'thickness' (strain maps)
%     sma_section_label(cfg,'suffix')  text appended to a title  ''  or  ' - SMA core, h = +-5 mm'
%     sma_section_label(cfg,'core')    'SMA core, h = +-5 mm' / 'h = +-5 mm'
if nargin < 2, kind = 'axis'; end
comp = isfield(cfg,'hle') && ~isempty(cfg.hle);
h = cfg.geom.thick/2;
core = sprintf('h = \\pm%g mm', h);
if comp, core = ['SMA core, ' core]; end
switch kind
    case 'axis',   s = 'h [mm]';  if comp, s = 'h [mm] (SMA core)'; end
    case 'thickness', s = 'thickness [mm]';  if comp, s = 'thickness [mm] (SMA core)'; end
    case 'suffix', s = ''; if comp, s = [' - ' core]; end
    case 'core',   s = core;
    otherwise, error('sma_section_label: unknown kind "%s"', kind);
end
end
