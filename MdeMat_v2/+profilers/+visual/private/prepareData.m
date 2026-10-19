function subdata = prepareData(dataset, varName, nDims)
% PREPAREDATA Prepares data for rendering from datasets.
%
% Syntax:
%   subdata = prepareData(dataset, varName, nDims)
%
% Inputs:
%   dataset - Map of dataset names to dataset objects
%   varName - Name of the variable to prepare
%   nDims - Number of spatial dimensions (2 or 3)
%
% Outputs:
%   subdata - Map of dataset names to prepared data structures

dataNames = fieldnames(dataset);
[fieldName, varIdx] = splitVarName(varName);
subdata = struct();

for i = 1:length(dataNames)
    dataName = dataNames{i};
    g = dataset.(dataName);

    if ~isfield(g.data, fieldName)
        continue;
    end

    if isempty(g.coords), continue; end

    x = g.coords;
    u = g.data.(fieldName)(:, varIdx);

    % Parse coordinates based on format and dimensions
    if iscell(x)
        coords = x(1:nDims);
    else
        coords = cell(1, nDims);
        for iDim = 1:nDims
            coords{iDim} = unique(x(iDim, :));
        end
    end

    % Reshape data if needed
    if nDims == 1
        n = [1, length(coords{1})];
    else
        n = cellfun(@length, coords);
    end
    if length(u) == prod(n)
        u = reshape(u, n(:).');
    end

    % Create data structure
    dataStruct = struct('u', u);
    for iDim = 1:nDims
        coordName = sprintf('x%d', iDim);
        dataStruct.(coordName) = coords{iDim};
    end

    subdata.(dataName) = dataStruct;
end
end