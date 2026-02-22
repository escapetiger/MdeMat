function [fieldName, varIdx] = splitVarName(str)
tokens = regexp(str, '^([a-zA-Z_]+)(\d+)$', 'tokens');
core.except.assert(~isempty(tokens), 'InvalidInput', ...
    'Input does not match the expected format.');
fieldName = tokens{1}{1};
varIdx = str2double(tokens{1}{2});
end