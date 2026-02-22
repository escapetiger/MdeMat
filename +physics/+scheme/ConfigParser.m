classdef ConfigParser < handle
    % CONFIGPARSER Configuration file parser for KEY=VALUE format.
    %
    %   ConfigParser provides a robust interface for reading and processing
    %   configuration files that use KEY=VALUE syntax. Supports extracting
    %   options as MATLAB variables and handling various configuration
    %   settings with default values and validation.
    %
    %   The parser supports configuration directives:
    %
    %   - OPTION=name=value: Defines configuration options with automatic
    %   type conversion
    %
    %   Type conversion is performed automatically for numeric values,
    %   logical values (true/false), and MATLAB expressions. String values
    %   are preserved as-is. Default values and validators can be
    %   configured for each option using the addConfig method.
    %
    %   Configuration files use KEY=VALUE format with # for comments. Empty
    %   lines and comment lines are ignored during parsing. All parsed
    %   options are stored in the Config property.
    %
    % See also:
    %   core.except.assert, core.except.verify

    properties
        Config % Struct containing all parsed configuration options
        Validators % Struct containing validator functions
    end

    methods
        function obj = ConfigParser()
            % CONFIGPARSER Construct an instance of ConfigParser.
            %
            %   obj = ConfigParser() creates a new ConfigParser instance
            %   with empty configuration and validators structures.

            obj.Config = struct();
            obj.Validators = struct();
        end

        function obj = loadFromFile(obj, file)
            % LOADFROMFILE Load configuration from file and add options to Config.
            %
            %   obj = loadFromFile(obj, file) reads and processes the
            %   specified configuration file @a file, extracting
            %   configuration options. All parsed options are added to obj.Config.

            arguments
                obj physics.scheme.ConfigParser
                file {mustBeTextScalar, mustBeFile}
            end

            %< Use static method to parse file
            fileOptions = physics.scheme.ConfigParser.parseFile(file);

            %< Transfer parsed options to Config with validation
            fieldNames = fieldnames(fileOptions);
            for i = 1:length(fieldNames)
                name = fieldNames{i};
                value = fileOptions.(name);
                obj.setConfig(name, value);
            end
        end

        function obj = addConfig(obj, name, options)
            % ADDCONFIG Add configuration option with default value and
            % validator.
            %
            %   obj = addConfig(obj, name, options) adds a configuration
            %   option with the specified @a name and @a options structure
            %   containing default value and validator function. The
            %   options structure should contain fields 'default' and
            %   'validator'.

            arguments
                obj physics.scheme.ConfigParser
                name {mustBeTextScalar}
                options.default = []
                options.validator = []
            end

            obj.Config.(name) = options.default;
            obj.Validators.(name) = options.validator;
        end

        function tf = hasConfig(obj, name)
            % HASCONFIG Check if configuration option exists.
            %
            %   tf = hasConfig(obj, name) returns true if the configuration
            %   option with the specified @a name exists in obj.Config,
            %   false otherwise.

            arguments
                obj physics.scheme.ConfigParser
                name {mustBeTextScalar}
            end

            tf = isfield(obj.Config, name);
        end

        function value = getConfig(obj, name)
            % GETCONFIG Get configuration option value.
            %
            %   value = getConfig(obj, name) returns the value of the
            %   configuration option with the specified @a name from
            %   obj.Config.

            arguments
                obj physics.scheme.ConfigParser
                name {mustBeTextScalar}
            end

            core.except.assert(isfield(obj.Config, name), 'ConfigNotFound', ...
                'Configuration option not found: %s', name);

            value = obj.Config.(name);
        end

        function obj = setConfig(obj, name, value)
            % SETCONFIG Set configuration option value.
            %
            %   obj = setConfig(obj, name, value) sets the configuration
            %   option with the specified @a name to the given @a value.
            %   If a validator exists for this option, it will be applied.

            arguments
                obj physics.scheme.ConfigParser
                name {mustBeTextScalar}
                value
            end

            if isfield(obj.Validators, name)
                validator = obj.Validators.(name);
                if isa(validator, 'function_handle')
                    core.except.assert(validator(value), 'InvalidConfig', ...
                        'Configuration parameter %s is invalid.', name);
                end
            end

            obj.Config.(name) = value;
        end
    end

    methods (Static)
        function config = parseFile(file)
            % PARSEFILE Parse configuration file and return extracted options.
            %
            %   config = parseFile(file) reads and processes the specified
            %   configuration file @a file, extracting configuration options.
            %
            %   The parser supports configuration directives:
            %
            %   - OPTION=name=value: Defines configuration options with automatic
            %   type conversion
            %
            %   Type conversion is performed automatically for numeric values,
            %   logical values (true/false), and MATLAB expressions. String values
            %   are preserved as-is.
            %
            %   Configuration files use KEY=VALUE format with # for comments. Empty
            %   lines and comment lines are ignored during parsing.

            arguments
                file {mustBeTextScalar, mustBeFile}
            end

            settings = physics.scheme.ConfigParser.parseFileImpl(file);
            config = struct();

            keys = fieldnames(settings);
            for i = 1:length(keys)
                key = keys{i};
                value = settings.(key);

                if strcmp(key, 'OPTION')
                    for j = 1:length(value)
                        config = physics.scheme.ConfigParser.addOption(value{j}, config);
                    end
                else
                    core.except.verify(0, 'UnknownSetting', ...
                        'Unknown setting: %s\n', key);
                end
            end
        end
    end

    methods (Static, Access = private)
        function settings = parseFileImpl(file)
            % PARSEFILEIMPL Parse configuration file into structured settings format.

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
                    core.except.verify(0, 'InvalidOption', ...
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
                'Invalid OPTION format (expected name=value): %s', ...
                optionName);

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
    end
end