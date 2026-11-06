function spec = primitiveSpec(name)
% The fixed v1 grammar has no class-name-driven construction.
names = {'double','single','logical','int8','uint8','int16','uint16','int32','uint32','int64','uint64'};
widths = [8,4,1,1,1,2,2,4,4,8,8];
if ~ischar(name) || ~isrow(name)
    error('mdemat:io:InvalidRecord','Primitive type must be a character row.');
end
index = find(strcmp(name,names),1);
if isempty(index), error('mdemat:io:InvalidRecord','Unsupported primitive type.'); end
spec = struct('Bytes',widths(index),'Complex',ismember(name,{'double','single'}), ...
    'Sparse',ismember(name,{'double','logical'}));
end
