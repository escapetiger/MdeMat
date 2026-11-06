function names = checkpointModelKeys(names)
%CHECKPOINTMODELKEYS Exclude known observation/output settings at execution boundary.
%   A model can instead pass explicit modelKeys to checkpointConfiguration.
excluded = {'timer','visualizer','analyzer','analysis','verbose','paths', ...
    'components','exact','errorTypes','eId','titlePrefix','title','plotOptions'};
names = names(~ismember(names,excluded) & ~startsWith(names,'ckpt'));
end
