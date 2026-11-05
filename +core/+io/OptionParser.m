classdef OptionParser < handle
    % OPTIONPARSER Configuration file parser for KEY=VALUE format.
    %
    %   OptionParser provides a robust interface for reading and processing
    %   configuration files that use KEY=VALUE syntax. Supports extracting
    %   options as MATLAB variables, dynamically adding directories to the
    %   MATLAB path, and handling various configuration settings.
    %
    % Examples:
    %   % Basic configuration parsing
    %   parser = core.io.OptionParser();
    %   options = parser.parse('config/simulation.conf');
    %
    %   % Parse with cleanup
    %   parser = core.io.OptionParser();
    %   options = parser.parse('config/setup.conf');
    %   parser.reset();  % Clean up paths
    %
    % See Also:
    %   core.except.assert, core.except.verify

    properties
        paths % Cell array of MATLAB paths added
    end

    methods
        function options = parse(obj, file)
            % PARSE Parse configuration file and return extracted options.
            %
            %   options = parse(obj, file) reads and processes the
            %   specified configuration file @a file, extracting
            %   configuration options and adding paths to MATLAB's search
            %   path as specified.
            %
            % Inputs:
            %   obj - The OptionParser object
            %   file - Path to the configuration file
            %
            % Outputs:
            %   options - Structure containing parsed configuration options

            core.except.assert(isfile(file), ...
                'FileNotFound', ...
                'Configuration file not found: %s', file);

            settings = obj.parseFile(file);
            options = struct();
            newPaths = {};

            keys = fieldnames(settings);
            for i = 1:length(keys)
                key = keys{i};
                value = settings.(key);

                if strcmp(key, 'OPTION')
                    for j = 1:length(value)
                        options = obj.addOption(value{j}, options);
                    end
                elseif strcmp(key, 'MATLABPATH')
                    newPaths = obj.addMatlabPath(value, newPaths);
                elseif strcmp(key, 'SUBDIRS')
                    if strcmp(value, ''), continue; end
                    newPaths = obj.addSubdirPath(value, newPaths);
                else
                    core.except.verify(0, 'UnknownSetting', ...
                        'Unknown setting: %s\n', key);
                end
            end

            if ~isempty(newPaths)
                addpath(strjoin(newPaths, pathsep));
                obj.paths = [obj.paths, newPaths];
            end
        end

        function obj = reset(obj)
            % RESET Remove all MATLAB paths added during parsing operations.
            %
            %   obj = reset(obj) removes all MATLAB search paths that were
            %   added by the parser during previous parse() operations.
            %
            % Inputs:
            %   obj - The OptionParser object
            %
            % Outputs:
            %   obj - The OptionParser object

            if ~isempty(obj.paths)
                rmpath(strjoin(obj.paths, pathsep));
                obj.paths = {};
                fprintf('Successfully reset option parser and removed %d paths.\n', ...
                    length(obj.paths));
            else
                fprintf('Option parser reset completed (no paths to remove).\n');
            end
        end
    end

    methods (Static, Access = private)
        function settings = parseFile(file)
            % PARSEFILE Parse configuration file into structured settings format.

            settings = struct();
            fid = fopen(file, 'r');
            core.except.assert(fid ~= -1, ...
                'FileOpenError', ...
                'Cannot open configuration file: %s', file);

            cleanupObj = onCleanup(@() fclose(fid));

            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if isempty(line) || startsWith(line, '#')
                    continue;
                end

                tokens = split(line, '=');
                if numel(tokens) < 2
                    core.except.verify(0, 'InvalidConfig', ...
                        'Invalid configuration line format: %s', line);
                    continue;
                end

                if strcmp(tokens{1}, 'OPTION') && numel(tokens) >= 3
                    key = 'OPTION';
                    valueTokens = tokens(2:end);
                    valueStr = strjoin(valueTokens, '=');
                    if isfield(settings, key)
                        settings.(key){end +1} = valueStr;
                    else
                        settings.(key) = {valueStr};
                    end
                else
                    key = strtrim(tokens{1});
                    valueTokens = tokens(2:end);
                    value = strjoin(valueTokens, '=');
                    settings.(key) = strtrim(value);
                end
            end
        end

        function options = addOption(optionName, options)
            % ADDOPTION Parse OPTION string and add to options structure.

            tokens = split(optionName, '=', 2);
            core.except.assert(numel(tokens) == 2, 'InvalidOption', ...
                'Invalid OPTION format (expected name=value): %s', optionName);

            name = strtrim(tokens{1});
            value0 = strtrim(tokens{2});

            value = str2double(value0);
            if isnan(value)
                if strcmpi(value0, 'true')
                    value = true;
                elseif strcmpi(value0, 'false')
                    value = false;
                else
                    try
                        value = eval(value0);
                    catch
                        value = value0;
                    end
                end
            end
            options.(name) = value;
        end

        function paths = addMatlabPath(pathString, paths)
            % ADDMATLABPATH Add directory path to paths collection.

            paths = [paths, {pathString}];
        end

        function paths = addSubdirPath(pathString, paths)
            % ADDSUBDIRPATH Recursively add subdirectories to paths collection.

            dirs = strsplit(pathString, ',');
            for i = 1:length(dirs)
                baseDir = strtrim(dirs{i});
                if isfolder(baseDir)
                    genPathStr = genpath(baseDir);
                    newPaths = strsplit(genPathStr, pathsep);
                    newPaths = newPaths(~cellfun('isempty', newPaths));
                    paths = [paths, newPaths];
                else
                    core.except.verify(0, 'DirectoryNotFound', ...
                        'Base directory not found: %s', baseDir);
                end
            end
        end
    end
end