function dataset = sampledDataset(data)
%SAMPLEDDATASET Adapt the retained sampled-study value layout at its boundary.
%   DATA.x is dimension-by-points; other fields are points-by-components.
%   A monotone 1D sample sequence is explicitly an ordered axis. Other point
%   sets remain points, without invented topology or independently sorted data.
arguments
    data (1,1) struct
end
fields=struct();
if isempty(fieldnames(data))
    dataset=struct('Type','fixed','Data',fields);return;
end
if ~isfield(data,'x')||~isnumeric(data.x)||~ismatrix(data.x)
    error('mdemat:execution:sampledDataset:Coordinates','Supply the saved physical points in x.');
end
points=data.x;names=setdiff(fieldnames(data),{'x'},'stable');
for index=1:numel(names)
    name=names{index};values=data.(name);
    if ~isnumeric(values)||~ismatrix(values)||size(values,1)~=size(points,2)
        error('mdemat:execution:sampledDataset:Values','Saved fields require one row per point.');
    end
    view=struct('Values',values.','ComponentIds',1:size(values,2));
    if size(points,1)==1 && (size(points,2)<2 || all(diff(points)>0)||all(diff(points)<0))
        view.Axes={points};
    else
        view.Coordinates=points;
    end
    fields.(name)=view;
end
dataset=struct('Type','fixed','Data',fields);
end
