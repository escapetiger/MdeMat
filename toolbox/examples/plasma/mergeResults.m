function mergeResults(inputDir, outputDir, groupBy, options)
% MERGERESULTS Merge CSV result files by specified parameter (vertical concatenation).
%
%   mergeResults(inputDir, outputDir, groupBy) scans @a inputDir for CSV
%   files matching the pattern LDI_FDGHS2_<metric>_epsilon<eps>_tau0<tau>.csv,
%   groups them by the specified parameter @a groupBy, and writes merged
%   tables to @a outputDir.
%
%   The merged tables vertically concatenate all data from files sharing the
%   same grouping parameter value, with additional columns for varying parameters.
%   Files within tau0 groups use descending epsilon; epsilon groups use
%   ascending tau0; metric groups use ascending [epsilon,tau0]. Source row
%   order and parameter spelling in output names are retained. Decimal and
%   scientific numeric parameter names are accepted. Existing merged files
%   are refused unless Overwrite=true is explicitly supplied.
%
% See also:
%   readtable, writetable

    arguments
        inputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')), 'ex01ldi')
        outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')), 'ex01ldi')
        groupBy (1,1) string {mustBeMember(groupBy, ["tau0", "epsilon", "metric"])} = "tau0"
        options.Overwrite (1,1) logical = false
    end

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    csvFiles = dir(fullfile(inputDir, '*.csv'));
    if isempty(csvFiles)
        warning('No CSV files found in %s', inputDir);
        return;
    end

    fileInfos = parseFileNames({csvFiles.name}');

    switch groupBy
        case "tau0"
            mergeByTau0(fileInfos, inputDir, outputDir, options.Overwrite);
        case "epsilon"
            mergeByEpsilon(fileInfos, inputDir, outputDir, options.Overwrite);
        case "metric"
            mergeByMetric(fileInfos, inputDir, outputDir, options.Overwrite);
    end

    fprintf('Merged results written to: %s\n', outputDir);
end

function fileInfos = parseFileNames(fileNames)
% PARSEFILENAMES Parse parameter values from file names.

    nFiles = numel(fileNames);
    fileInfos = struct('fileName', cell(nFiles, 1), ...
                       'metric', cell(nFiles, 1), ...
                       'epsilon', cell(nFiles, 1), ...
                       'tau0', cell(nFiles, 1), ...
                       'epsilonStr', cell(nFiles, 1), ...
                       'tau0Str', cell(nFiles, 1));

    number = '[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?';
    pattern = ['^LDI_FDGHS2_(\w+)_epsilon(', number, ')_tau0(', number, ')\.csv$'];

    for i = 1:nFiles
        tokens = regexp(fileNames{i}, pattern, 'tokens', 'once');
        if ~isempty(tokens)
            fileInfos(i).fileName = fileNames{i};
            fileInfos(i).metric = tokens{1};
            fileInfos(i).epsilonStr = tokens{2};
            fileInfos(i).tau0Str = tokens{3};
            fileInfos(i).epsilon = str2double(tokens{2});
            fileInfos(i).tau0 = str2double(tokens{3});
        end
    end

    validIdx = ~cellfun(@isempty, {fileInfos.fileName});
    fileInfos = fileInfos(validIdx);
end

function mergeByTau0(fileInfos, inputDir, outputDir, overwrite)
% MERGEBYTAU0 Merge files with the same tau0 value vertically.

    uniqueTau0 = unique([fileInfos.tau0]);
    uniqueMetrics = unique({fileInfos.metric});

    for iTau = 1:numel(uniqueTau0)
        tau0Val = uniqueTau0(iTau);
        tau0Str = fileInfos(find([fileInfos.tau0] == tau0Val, 1)).tau0Str;

        for iMetric = 1:numel(uniqueMetrics)
            metricName = uniqueMetrics{iMetric};

            idx = ([fileInfos.tau0] == tau0Val) & ...
                  strcmp({fileInfos.metric}, metricName);
            subset = fileInfos(idx);

            if isempty(subset)
                continue;
            end

            [~, sortIdx] = sort([subset.epsilon], 'descend');
            subset = subset(sortIdx);

            allTables = cell(numel(subset), 1);
            for iFile = 1:numel(subset)
                filePath = fullfile(inputDir, subset(iFile).fileName);
                data = readtable(filePath);
                nRows = height(data);
                data.epsilon = repmat(subset(iFile).epsilon, nRows, 1);
                data = movevars(data, 'epsilon', 'Before', 1);
                allTables{iFile} = data;
            end

            mergedTable = vertcat(allTables{:});

            outFileName = sprintf('LDI_FDGHS2_%s_tau0%s_merged.csv', ...
                                  metricName, tau0Str);
            writeMergedTable(mergedTable, fullfile(outputDir, outFileName), overwrite);
            fprintf('  Created: %s\n', outFileName);
        end
    end
end

function mergeByEpsilon(fileInfos, inputDir, outputDir, overwrite)
% MERGEBYEPSILON Merge files with the same epsilon value vertically.

    uniqueEps = unique([fileInfos.epsilon]);
    uniqueMetrics = unique({fileInfos.metric});

    for iEps = 1:numel(uniqueEps)
        epsVal = uniqueEps(iEps);
        epsStr = fileInfos(find([fileInfos.epsilon] == epsVal, 1)).epsilonStr;

        for iMetric = 1:numel(uniqueMetrics)
            metricName = uniqueMetrics{iMetric};

            idx = ([fileInfos.epsilon] == epsVal) & ...
                  strcmp({fileInfos.metric}, metricName);
            subset = fileInfos(idx);

            if isempty(subset)
                continue;
            end

            [~, sortIdx] = sort([subset.tau0]);
            subset = subset(sortIdx);

            allTables = cell(numel(subset), 1);
            for iFile = 1:numel(subset)
                filePath = fullfile(inputDir, subset(iFile).fileName);
                data = readtable(filePath);
                nRows = height(data);
                data.tau0 = repmat(subset(iFile).tau0, nRows, 1);
                data = movevars(data, 'tau0', 'Before', 1);
                allTables{iFile} = data;
            end

            mergedTable = vertcat(allTables{:});

            outFileName = sprintf('LDI_FDGHS2_%s_epsilon%s_merged.csv', ...
                                  metricName, epsStr);
            writeMergedTable(mergedTable, fullfile(outputDir, outFileName), overwrite);
            fprintf('  Created: %s\n', outFileName);
        end
    end
end

function mergeByMetric(fileInfos, inputDir, outputDir, overwrite)
% MERGEBYMETRIC Merge files with the same metric type vertically.

    uniqueMetrics = unique({fileInfos.metric});

    for iMetric = 1:numel(uniqueMetrics)
        metricName = uniqueMetrics{iMetric};

        idx = strcmp({fileInfos.metric}, metricName);
        subset = fileInfos(idx);

        if isempty(subset)
            continue;
        end

        [~, sortIdx] = sortrows([[subset.epsilon]', [subset.tau0]']);
        subset = subset(sortIdx);

        allTables = cell(numel(subset), 1);
        for iFile = 1:numel(subset)
            filePath = fullfile(inputDir, subset(iFile).fileName);
            data = readtable(filePath);
            nRows = height(data);
            data.epsilon = repmat(subset(iFile).epsilon, nRows, 1);
            data.tau0 = repmat(subset(iFile).tau0, nRows, 1);
            data = movevars(data, {'epsilon', 'tau0'}, 'Before', 1);
            allTables{iFile} = data;
        end

        mergedTable = vertcat(allTables{:});

        outFileName = sprintf('LDI_FDGHS2_%s_all_merged.csv', metricName);
        writeMergedTable(mergedTable, fullfile(outputDir, outFileName), overwrite);
        fprintf('  Created: %s\n', outFileName);
    end
end

function writeMergedTable(data, file, overwrite)
    assert(overwrite || ~isfile(file), 'mdemat:examples:FileExists', ...
        'Merged output already exists: %s. Set Overwrite=true to replace it.', file);
    writetable(data, file);
end
