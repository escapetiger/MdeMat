function outId = backtrack(id)
% BACKTRACK Construct a identifier with module context.
%
%   outId = backtrack(id) extracts context information from the call stack
%   to build hierarchical identifiers of form 'module:identifier'.
%
% Inputs:
%   id - Base identifier string (e.g., 'InvalidInput')
%
% Outputs:
%   outId - Full identifier with context (e.g., 'core:except:InvalidInput')
%
% Examples:
%   % Basic usage - construct contextual identifier
%   id = backtrack('InvalidInput');
%   % Returns something like 'core:except:InvalidInput'
%
%   % Usage in error handling
%   fullId = backtrack('OutOfRange');
%   error(fullId, 'Parameter value is out of acceptable range.');
%
% See Also:
%   core.except.assert, core.except.verify

LIBRARY = 'MdeMat';

stack = dbstack(2, '-completenames');
if isempty(stack)
    caller = 'base';
    module = 'unknown';
else
    [~, caller, ~] = fileparts(stack(1).name);
    
    path = strsplit(stack(1).file, filesep);
    module = 'unknown';
    
    root = find(strcmp(path, LIBRARY), 1);
    
    if ~isempty(root) && root < length(path)
        module = '';
        for i = root + 1:length(path)
            if startsWith(path{i}, '+')
                if ~isempty(module)
                    module = strjoin({module, path{i}(2:end)}, ':');
                else
                    module = path{i}(2:end);
                end
            end
        end
    end
end

outId = [module, ':', caller, ':', id];

end