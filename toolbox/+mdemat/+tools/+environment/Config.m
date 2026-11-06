classdef Config
    %CONFIG Value map from declared keys to values and validation predicates.
    %   cfg = cfg.define(key, value, validator) declares a complete entry.
    %   Definitions may contain deferred defaults. set/update validate written
    %   values; validate(keys) checks the keys needed by the current operation.
    %   A validator returns one logical scalar and must not mutate external
    %   state. [] disables validation. Values containing handles remain shared,
    %   just as they do inside a MATLAB struct; no recursive copy is attempted.
    %
    %   All updates return a new Config and require explicit assignment. File
    %   parsing, cross-key model constraints and persistence belong to callers.
    %   See also mdemat.tools.environment.readConfig.

    properties (Access = private)
        Entries (1, 1) struct = struct()
    end

    methods
        function obj = define(obj, key, value, validator, options)
            %DEFINE Declare a value/predicate pair; override must be explicit.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                key {mustBeTextScalar}
                value
                validator = []
                options.override (1, 1) logical = false
            end
            key = configKey(key);
            if isfield(obj.Entries, key) && ~options.override
                error('mdemat:tools:environment:Config:DuplicateKey', ...
                    'Configuration key "%s" is already defined.', key);
            end
            if ~isequal(validator, []) && ~isa(validator, 'function_handle')
                error('mdemat:tools:environment:Config:InvalidValidator', ...
                    'The validator for "%s" must be a predicate or [].', key);
            end
            entry.Value = value;
            entry.Validator = validator;
            obj.Entries.(key) = entry;
        end

        function [value, validator] = get(obj, key)
            %GET Return a stored value and, optionally, its predicate.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                key {mustBeTextScalar}
            end
            key = configKey(key);
            obj.requireKey(key);
            value = obj.Entries.(key).Value;
            validator = obj.Entries.(key).Validator;
        end

        function tf = has(obj, key)
            %HAS Return a logical scalar indicating whether key is declared.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                key {mustBeTextScalar}
            end
            tf = isfield(obj.Entries, configKey(key));
        end

        function names = keys(obj)
            %KEYS Return declared names as a column cell array of char vectors.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
            end
            names = fieldnames(obj.Entries);
        end

        function obj = set(obj, key, value)
            %SET Replace one value without changing its predicate.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                key {mustBeTextScalar}
                value
            end
            key = configKey(key);
            obj.requireKey(key);
            obj.checkValue(key, value);
            obj.Entries.(key).Value = value;
        end

        function obj = update(obj, values)
            %UPDATE Validate a scalar struct patch before replacing any values.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                values (1, 1) struct
            end
            names = fieldnames(values);
            for index = 1:numel(names)
                obj.requireKey(names{index});
                obj.checkValue(names{index}, values.(names{index}));
            end
            for index = 1:numel(names)
                obj.Entries.(names{index}).Value = values.(names{index});
            end
        end

        function validate(obj, names)
            %VALIDATE Check selected keys; omitted names checks every entry.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
                names = []
            end
            if nargin < 2, names = obj.keys(); end
            names = configKeys(names);
            for index = 1:numel(names)
                obj.requireKey(names{index});
                obj.checkValue(names{index}, obj.Entries.(names{index}).Value);
            end
        end

        function values = toStruct(obj)
            %TOSTRUCT Return values only, retaining MATLAB's shallow semantics.
            arguments
                obj (1, 1) mdemat.tools.environment.Config
            end
            values = struct();
            names = obj.keys();
            for index = 1:numel(names)
                values.(names{index}) = obj.Entries.(names{index}).Value;
            end
        end
    end

    methods (Access = private)
        function requireKey(obj, key)
            if ~isfield(obj.Entries, key)
                error('mdemat:tools:environment:Config:UnknownKey', ...
                    'Configuration key "%s" has not been declared.', key);
            end
        end

        function checkValue(obj, key, value)
            predicate = obj.Entries.(key).Validator;
            if isempty(predicate), return; end
            try
                valid = predicate(value);
            catch exception
                failure = MException('mdemat:tools:environment:Config:InvalidValue', ...
                    'The validator for configuration key "%s" failed.', key);
                throwAsCaller(addCause(failure, exception));
            end
            if ~islogical(valid) || ~isscalar(valid)
                error('mdemat:tools:environment:Config:InvalidValidatorResult', ...
                    'The validator for "%s" must return one logical scalar.', key);
            end
            if ~valid
                error('mdemat:tools:environment:Config:InvalidValue', ...
                    'Configuration value "%s" is invalid.', key);
            end
        end
    end
end

function key = configKey(key)
key = char(key);
if ~isvarname(key)
    error('mdemat:tools:environment:Config:InvalidKey', ...
        'Configuration keys must be valid MATLAB identifiers.');
end
end

function names = configKeys(names)
if isempty(names) && ~ischar(names), names = {}; return; end
if ischar(names)
    names = {names};
elseif isstring(names)
    names = cellstr(names(:));
end
if ~iscell(names) || ~all(cellfun(@(name) ...
        (ischar(name) && isrow(name)) || (isstring(name) && isscalar(name)), names), 'all')
    error('mdemat:tools:environment:Config:InvalidKey', ...
        'Key selections must contain text scalars.');
end
names = cellfun(@configKey, names(:), 'UniformOutput', false);
end
