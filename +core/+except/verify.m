function verify(flag, id, msg, varargin)
% VERIFY Advanced verification.
%
%   verify(flag) throws a warning if @a flag is false. 
%
%   verify(flag, id, msg, ...) throws a warning with the specified @a id
%   and @a msg if @a flag is false.
%
% Inputs:
%   flag     - Logical condition to check
%   id       - Warning identifier (e.g., 'InvalidInput' or 'OutOfRange')
%   msg      - Warning message format string
%   varargin - Optional variables for the format string
%
% Outputs:
%   NULL
%
% Examples:
%   % Basic usage with simple condition
%   core.except.verify(x > 0, 'NonPositive', 'Value should be positive.');
%
%   % With formatted message
%   core.except.verify(n < 100, 'TooLarge', 'Value %d exceeds limit.', n);
%
%   % Default behavior (minimal arguments)
%   core.except.verify(isreal(x));
%
% See Also:
%   core.except.assert, core.except.backtrack

if nargin == 1
    if ~flag, warning('Verification failed.'); end
    return;
end

if ~flag
    outId = backtrack(id);
    
    if ~isempty(varargin)
        outMsg = sprintf(msg, varargin{:});
    else
        outMsg = msg;
    end
    
    warning(outId, outMsg);
end

end