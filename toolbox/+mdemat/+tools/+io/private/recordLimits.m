function limits = recordLimits(overrides)
% Resource policy shared by readers and writers; all limits are adjustable.
limits = struct('MaxFileBytes',2^31,'MaxArrayBytes',2^31-1,'MaxTotalBytes',2^32, ...
    'MaxElements',1e8,'MaxNodes',1e6,'MaxDepth',128,'MaxDimensions',32);
if nargin == 0, return; end
if ~isstruct(overrides) || ~isscalar(overrides) || ...
        ~all(ismember(fieldnames(overrides),fieldnames(limits)))
    error('mdemat:io:InvalidLimits','Use a scalar struct of documented resource limits.');
end
names = fieldnames(overrides);
for index = 1:numel(names)
    value = overrides.(names{index});
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value < 1 || value ~= fix(value)
        error('mdemat:io:InvalidLimits','Resource limits must be finite positive integers.');
    end
    limits.(names{index}) = double(value);
end
end
