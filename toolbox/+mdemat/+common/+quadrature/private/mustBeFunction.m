function mustBeFunction(x)
% MUSTBEFUNCTION Validation function for function handles and Function objects.
%
%   mustBeFunction(x) validates that x is either a function handle or a
%   mdemat.data.function.Function object. Throws an error if the validation fails.

if ~isa(x, 'mdemat.data.function.Function') && ~isa(x, 'function_handle')
    error('Must be mdemat.data.function.Function or function handle.');
end
end