function verify(flag, id, msg, varargin)
% VERIFY Advanced verification with contextual warning reporting.
%
%   verify(flag) throws a warning if @a flag is false.
%
%   verify(flag, id, msg, ...) throws a warning with the specified @a id
%   and @a msg if @a flag is false. The @a id is automatically prefixed with
%   module context information. Additional arguments are used for string
%   formatting of @a msg using sprintf.
%
% See also:
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