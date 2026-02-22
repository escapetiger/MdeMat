function mustBeFunction(x)
% MUSTBEFUNCTION Validation function for function handles and Function objects.
%
%   mustBeFunction(x) validates that x is either a function handle or a
%   core.function.Function object. Throws an error if the validation fails.

if ~isa(x, 'core.function.Function') && ~isa(x, 'function_handle')
    error('Must be core.function.Function or function handle.');
end
end