function S = cut(T, d)
% CUT Cuts a tensor along a specified dimension.
%
%   S = cut(T, d) returns slices of @a T along the dimension @a d.
%
% Inputs:
%   T - Input tensor
%   d - Dimension index
%
% Outputs:
%   S - Cell array of slices of @a T along the dimension @a d
%
% Examples:
%   % Create a 3D T and cut along the second dimension
%   T = rand(3, 4, 5);
%   S = cut(T, 2);
%   % returns a cell array of length 4 of which entry is a 3×5 matrix

if isempty(T)
    S = {};
    return;
end
p = [d, setdiff(1:ndims(T), d)];
T = permute(T, p);
S = num2cell(T, 2:ndims(T));
S = cellfun(@(s) squeeze(ipermute(s, p)), S, 'Un', 0);
end