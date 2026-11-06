function configurePlotLayouts(scheme, datasets)
%CONFIGUREPLOTLAYOUTS Bind stable field/component/dataset identities to figures.
%   Sampling and layout selection belong to execution. The graphics session
%   receives only value specifications and does not infer physical fields.
options = scheme.Observation.Plot;
if isempty(options.Modes), return; end
session = scheme.getObserver('visualizer');
if isempty(session)
    session = mdemat.tools.visual.Visualizer();
    scheme.setObserver('visualizer', session);
end
if ~isempty(session.Layouts), return; end
ids = string(fieldnames(options.DatasetTypes)).';
if isempty(ids), ids = string(fieldnames(datasets)).'; end
fields = fieldnames(options.Components);
layouts = struct('Name', {}, 'Size', {}, 'Panels', {});
modes = string(options.Modes);
for modeIndex = 1:numel(modes)
    mode = char(modes(modeIndex));
    panels = struct('Mode', {}, 'Field', {}, 'ComponentId', {}, 'DatasetIds', {});
    for field = 1:numel(fields)
        name = fields{field};
        components = options.Components.(name);
        for component = components(:).'
            panels(end+1) = struct('Mode', mode, 'Field', name, ...
                'ComponentId', component, 'DatasetIds', ids); %#ok<AGROW>
        end
    end
    if isempty(panels), continue; end
    if ismember(mode, {'2d', '3d'}), columns = max(1, numel(ids)); else, columns = 1; end
    layouts(end+1) = struct('Name', sprintf('profile%d', modeIndex), ...
        'Size', [numel(panels), columns], 'Panels', panels); %#ok<AGROW>
end
session.setLayouts(layouts);
end
