function S = cut(T, d)
% CUT Cuts a tensor along a specified dimension.
%
%   S = cut(T, d) returns slices of @a T along the dimension @a d.
%
%   This function cuts tensor @a T along dimension @a d, returning a cell
%   array where each entry contains a slice of the original tensor. The
%   operation preserves all dimensions except @a d, which is sliced.
%   Returns empty cell array if @a T is empty.

if isempty(T)
    S = {};
    return;
end
p = [d, setdiff(1:ndims(T), d)];
T = permute(T, p);
S = num2cell(T, 2:ndims(T));
S = cellfun(@(s) squeeze(ipermute(s, p)), S, 'Un', 0);
end