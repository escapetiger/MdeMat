function mustBeFunctionHandleOrEmpty(x)
% MUSTBEFUNCTIONHANDLEOREMPTY Validation function for function handles.
%
%   mustBeFunctionHandleOrEmpty(x) validates that x is either a function
%   handle or empty. Throws an error if the validation fails.

if ~isempty(x) && ~isa(x, 'function_handle')
    error('Must be empty or function handle.');
end
end