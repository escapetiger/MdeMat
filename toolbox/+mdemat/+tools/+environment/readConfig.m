function [values, info] = readConfig(file)
%READCONFIG Read explicit values without executing configuration code.
%   [values, info] = readConfig(file) parses OPTION=name=literal, MATLABPATH
%   and SUBDIRS. Paths are source-relative descriptions in info.Paths; path
%   and cwd are never changed. info.Lines maps OPTION names to source lines.
%   Duplicate options are rejected. Supported literals are quoted text,
%   logical/numeric scalars, rectangular arrays/cells and allowlisted scalar
%   arithmetic. Use parentheses around arithmetic entries containing spaces.
%   See also mdemat.tools.environment.Config.
arguments
    file {mustBeTextScalar, mustBeFile}
end
[ok, attributes] = fileattrib(char(file));
if ~ok
    error('mdemat:tools:environment:readConfig:FileOpenError', ...
        'Cannot resolve configuration file: %s', file);
end
file = attributes.Name;
[stream, message] = fopen(file, 'rt', 'n', 'UTF-8');
if stream < 0
    error('mdemat:tools:environment:readConfig:FileOpenError', '%s', message);
end
guard = onCleanup(@() fclose(stream));
values = struct();
info = struct('File', file, 'Lines', struct(), 'Paths', {{}});
lineNumber = 0;
while true
    line = fgetl(stream);
    if ~ischar(line)
        [message, code] = ferror(stream);
        if code ~= 0 && ~feof(stream)
            error('mdemat:tools:environment:readConfig:ReadFailed', ...
                'Cannot read %s: %s', file, message);
        end
        break;
    end
    lineNumber = lineNumber + 1;
    if lineNumber == 1 && ~isempty(line) && line(1) == char(65279)
        line(1) = [];
    end
    try
        line = strtrim(stripComment(line));
        if isempty(line), continue; end
        [directive, text] = splitAssignment(line);
        switch upper(directive)
            case 'OPTION'
                [name, literal] = splitAssignment(text);
                name = validName(name);
                if isfield(values, name)
                    error('mdemat:tools:environment:readConfig:DuplicateOption', ...
                        'Configuration option "%s" is repeated.', name);
                end
                values.(name) = parseLiteral(literal);
                info.Lines.(name) = lineNumber;
            case {'MATLABPATH', 'SUBDIRS'}
                directory = strtrim(text);
                if ~isempty(directory) && any(directory(1) == ['''', '"'])
                    directory = parseLiteral(directory);
                end
                if ~((ischar(directory) && isrow(directory)) || ...
                        (isstring(directory) && isscalar(directory)))
                    error('mdemat:tools:environment:readConfig:InvalidDirectory', ...
                        'Resource directories must be nonempty text scalars.');
                end
                directory = char(directory);
                if isempty(directory)
                    error('mdemat:tools:environment:readConfig:InvalidDirectory', ...
                        'Resource directories must be nonempty text scalars.');
                end
                if ~isAbsoluteDirectory(directory)
                    directory = fullfile(fileparts(file), directory);
                end
                if ~isfolder(directory)
                    error('mdemat:tools:environment:readConfig:DirectoryNotFound', ...
                        'Configuration directory does not exist: %s', directory);
                end
                [~, attributes] = fileattrib(directory);
                directory = attributes.Name;
                if strcmpi(directive, 'SUBDIRS')
                    directories = strsplit(genpath(directory), pathsep);
                    directories = directories(~cellfun('isempty', directories));
                else
                    directories = {directory};
                end
                info.Paths = unique([info.Paths, directories], 'stable');
            otherwise
                error('mdemat:tools:environment:readConfig:UnknownSetting', ...
                    'Unknown configuration directive: %s', directive);
        end
    catch exception
        identifier = exception.identifier;
        if ~startsWith(identifier, 'mdemat:tools:environment:readConfig:')
            identifier = 'mdemat:tools:environment:readConfig:InvalidConfiguration';
        end
        failure = MException(identifier, '%s:%d: %s', file, lineNumber, exception.message);
        throwAsCaller(addCause(failure, exception));
    end
end
end

function tf = isAbsoluteDirectory(directory)
if ispc
    tf = ~isempty(regexp(directory, '^[A-Za-z]:[\\/]|^[\\/]', 'once'));
else
    tf = startsWith(directory, filesep);
end
end

function name = validName(name)
name = char(name);
if ~isvarname(name)
    error('mdemat:tools:environment:readConfig:InvalidOption', ...
        'Invalid configuration option name: %s', name);
end
end

function [name, value] = splitAssignment(text)
separator = find(text == '=', 1);
if isempty(separator) || separator == 1
    error('mdemat:tools:environment:readConfig:InvalidOption', ...
        'Expected name=value, received: %s', text);
end
name = strtrim(text(1:separator-1));
value = strtrim(text(separator+1:end));
end

function text = stripComment(text)
quote = char(0);
index = 1;
while index <= length(text)
    symbol = text(index);
    if quote ~= char(0)
        if symbol == quote
            if index < length(text) && text(index+1) == quote
                index = index + 1;
            else
                quote = char(0);
            end
        end
    elseif any(symbol == ['''', '"'])
        quote = symbol;
    elseif symbol == '#'
        text = text(1:index-1);
        return;
    end
    index = index + 1;
end
end

function value = parseLiteral(text)
text = strtrim(text);
if isempty(text), value = ''; return; end
if any(text(1) == ['''', '"'])
    quote = text(1);
    if length(text) < 2 || text(end) ~= quote
        unsupported(text);
    end
    inner = text(2:end-1);
    % A quote inside the payload must be escaped by a doubled quote.
    bare = strrep(inner, [quote, quote], '');
    if contains(bare, quote), unsupported(text); end
    value = strrep(inner, [quote, quote], quote);
    if quote == '"'
        value = string(value);
    elseif isempty(value)
        value = '';
    end
elseif text(1) == '[' || text(1) == '{'
    isCell = text(1) == '{';
    if (isCell && text(end) ~= '}') || (~isCell && text(end) ~= ']')
        unsupported(text);
    end
    if isempty(strtrim(text(2:end-1)))
        if isCell, value = {}; else, value = []; end
        return;
    end
    rows = splitTopLevel(text(2:end-1), ';', false);
    blocks = cell(numel(rows), 1);
    for i = 1:numel(rows)
        columns = splitTopLevel(rows{i}, ',', true);
        entries = cellfun(@parseLiteral, columns, 'UniformOutput', false);
        if isCell
            blocks{i} = entries;
        else
            if ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && isscalar(x), entries))
                unsupported(text);
            end
            blocks{i} = [entries{:}];
        end
    end
    if numel(unique(cellfun(@(x) size(x, 2), blocks))) ~= 1
        unsupported(text);
    end
    value = vertcat(blocks{:});
elseif strcmpi(text, 'true')
    value = true;
elseif strcmpi(text, 'false')
    value = false;
elseif ~isempty(regexp(text, '^[+\-\d.]', 'once')) || ...
        ~isempty(regexp(text, '^(pi|eps|inf|nan)(\W|$)', 'once', 'ignorecase')) || ...
        contains(text, '(')
    value = parseArithmetic(text);
else
    if ~isempty(regexp(text, '[\[\]{};=]', 'once')), unsupported(text); end
    value = text;
end
end

function values = splitTopLevel(text, separator, whitespace)
values = {};
brackets = '';
quote = char(0);
start = 1;
index = 1;
afterSeparator = false;
while index <= length(text)
    symbol = text(index);
    if quote ~= char(0)
        if symbol == quote
            if index < length(text) && text(index+1) == quote
                index = index + 1;
            else
                quote = char(0);
            end
        end
    elseif any(symbol == ['''', '"'])
        quote = symbol;
    elseif any(symbol == '[{(')
        brackets(end+1) = symbol; %#ok<AGROW>
    elseif any(symbol == ']})')
        opening = '[{(';
        closing = ']})';
        if isempty(brackets) || brackets(end) ~= opening(closing == symbol)
            unsupported(text);
        end
        brackets(end) = [];
    elseif isempty(brackets) && (symbol == separator || (whitespace && isspace(symbol)))
        entry = strtrim(text(start:index-1));
        if ~isempty(entry)
            values{end+1} = entry; %#ok<AGROW>
            afterSeparator = false;
        end
        if symbol == separator
            if (isempty(entry) && (afterSeparator || isempty(values)))
                unsupported(text);
            end
            afterSeparator = true;
        end
        start = index + 1;
    end
    index = index + 1;
end
if quote ~= char(0) || ~isempty(brackets), unsupported(text); end
entry = strtrim(text(start:end));
if ~isempty(entry)
    values{end+1} = entry;
elseif afterSeparator || isempty(values)
    unsupported(text);
end
end

function value = parseArithmetic(input)
% Recursive descent over scalar constants and an explicit pure-function set.
text = input;
text = strrep(strrep(strrep(text, '.*', '*'), './', '/'), '.^', '^');
position = 1;
value = sumTerm();
skipSpaces();
if position <= length(text), unsupported(input); end
    function value = sumTerm()
        value = product();
        skipSpaces();
        while position <= length(text) && any(text(position) == '+-')
            op = text(position); position = position + 1;
            rhs = product();
            if op == '+', value = value + rhs; else, value = value - rhs; end
            skipSpaces();
        end
    end
    function value = product()
        value = unary();
        skipSpaces();
        while position <= length(text) && any(text(position) == '*/')
            op = text(position); position = position + 1;
            rhs = unary();
            if op == '*', value = value * rhs; else, value = value / rhs; end
            skipSpaces();
        end
    end
    function value = unary()
        skipSpaces();
        if position <= length(text) && any(text(position) == '+-')
            negative = text(position) == '-'; position = position + 1;
            value = unary();
            if negative, value = -value; end
        else
            value = powerTerm();
        end
    end
    function value = powerTerm()
        value = primary();
        skipSpaces();
        while position <= length(text) && text(position) == '^'
            position = position + 1;
            value = value ^ signedPrimary();
            skipSpaces();
        end
    end
    function value = signedPrimary()
        skipSpaces();
        if position <= length(text) && any(text(position) == '+-')
            negative = text(position) == '-'; position = position + 1;
            value = signedPrimary();
            if negative, value = -value; end
        else
            value = primary();
        end
    end
    function value = primary()
        skipSpaces();
        if position > length(text), unsupported(input); end
        if text(position) == '('
            position = position + 1; value = sumTerm();
            expect(')'); return;
        end
        number = regexp(text(position:end), ...
            '^(\d+(\.\d*)?|\.\d+)([eE][+\-]?\d+)?', 'match', 'once');
        if ~isempty(number)
            value = str2double(number); position = position + length(number); return;
        end
        name = regexp(text(position:end), '^[A-Za-z]+', 'match', 'once');
        if isempty(name), unsupported(input); end
        position = position + length(name);
        constants = struct('pi', pi, 'eps', eps, 'inf', Inf, 'nan', NaN);
        if isfield(constants, lower(name))
            value = constants.(lower(name)); return;
        end
        functions = struct('sin', @sin, 'cos', @cos, 'tan', @tan, ...
            'sqrt', @sqrt, 'exp', @exp, 'log', @log, 'abs', @abs);
        if ~isfield(functions, name), unsupported(input); end
        expect('('); argument = sumTerm(); expect(')');
        operation = functions.(name);
        value = operation(argument);
    end
    function expect(symbol)
        skipSpaces();
        if position > length(text) || text(position) ~= symbol, unsupported(input); end
        position = position + 1;
    end
    function skipSpaces()
        while position <= length(text) && isspace(text(position))
            position = position + 1;
        end
    end
end

function unsupported(text)
error('mdemat:tools:environment:readConfig:UnsupportedExpression', ...
    'Unsupported configuration literal/expression: %s', text);
end

