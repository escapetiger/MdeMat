function configureObservation(scheme)
%CONFIGUREOBSERVATION Supply model defaults without creating numerical tools.
%   Called once after candidate normalization. Explicit observation records
%   supplied by the caller are preserved; initialize may provide additional
%   model-owned physical reconstruction callbacks.
config = scheme.Config;
components = struct(); exacts = struct();
if isfield(config,'components'),components=config.components;end
if isfield(config,'exact') && ~isempty(config.exact)
    if isstruct(config.exact)
        exacts=config.exact;
    elseif isfield(components,'F')
        exacts.F=config.exact;
    else
        exacts.U=config.exact;
    end
end
plot=scheme.Observation.Plot; analysis=scheme.Observation.Analysis;
if isempty(fieldnames(plot.Components)),plot.Components=components;end
if isempty(fieldnames(analysis.Components)),analysis.Components=components;end
if isempty(fieldnames(plot.Exacts)),plot.Exacts=exacts;end
if isempty(fieldnames(analysis.Exacts)),analysis.Exacts=exacts;end
if isempty(plot.Title) && isfield(config,'titlePrefix'),plot.Title=config.titlePrefix;end
if isempty(fieldnames(analysis.Reductions))
    norms={'L1','L2','max'};
    if isfield(config,'errorTypes'),norms=config.errorTypes;end
    fields=fieldnames(analysis.Components);
    for index=1:numel(fields)
        analysis.Reductions.(fields{index})=struct('x',{norms},'v','');
    end
end
scheme.setObservation('plot',plot);scheme.setObservation('analysis',analysis);
end
