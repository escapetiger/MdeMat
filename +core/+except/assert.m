function assert(flag, id, msg, varargin)
% ASSERT Advanced assertion with contextual error reporting.
%
%   assert(flag) throws an error if @a flag is false.
%
%   assert(flag, id, msg, ...) throws an error with the specified @a id and
%   @a msg if @a flag is false. The @a id is automatically prefixed with
%   module context information. Additional arguments are used for string
%   formatting of @a msg using sprintf.
%
% See also:
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