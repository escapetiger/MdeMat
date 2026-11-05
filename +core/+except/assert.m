function assert(flag, id, msg, varargin)
% ASSERT Advanced assertion.
%
%   assert(flag) throws an error if @a flag is false. 
%
%   assert(flag, id, msg, ...) throws an error with the specified @a id and
%   @a msg if @a flag is false.
%
% Inputs:
%   flag     - Logical condition to check
%   id       - Error identifier (e.g., 'InvalidInput' or 'OutOfRange')
%   msg      - Error message format string
%   varargin - Optional variables for the format string
%
% Outputs:
%   NULL
%
% Examples:
%   % Basic usage with simple condition
%   core.except.assert(x > 0, 'NonPositive', 'Value must be positive.');
%
%   % With formatted message
%   core.except.assert(n < 100, 'TooLarge', 'Value %d exceeds limit.', n);
%
%   % Default behavior (minimal arguments)
%   core.except.assert(isreal(x));
%
% See Also:
%   core.except.verify, core.except.backtrack

if nargin == 1
    if ~flag, error('Assertion failed.'); end
    return;
end

if ~flag
    outId = backtrack(id);
    
    if ~isempty(varargin)
        outMsg = sprintf(msg, varargin{:});
    else
        outMsg = msg;
    end
    
    error(outId, outMsg);
end

end