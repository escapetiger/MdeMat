function varNames = parseVarNames(dataset)
% PARSEVARNAMES Collects all variable names from datasets.
%
% Syntax:
%   varNames = parseVarNames(dataset)
%
% Inputs:
%   dataset - Map of dataset names to dataset objects
%
% Outputs:
%   varNames - Cell array of all unique variable names

dataNames = fieldnames(dataset);
varNames = {};

for i = 1:length(dataNames)
    dataName = dataNames{i};
    data = dataset.(dataName);
    subVarNames = fieldnames(data);

    for j = 1:length(subVarNames)
        if ~ismember(subVarNames{j}, varNames)
            varNames{end+1} = subVarNames{j};
        end
    end
end
end