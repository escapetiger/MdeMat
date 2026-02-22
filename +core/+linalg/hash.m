function h = hash(A, opt)
%HASH Compute a hash value for a matrix (dense or sparse).
% 
%   H = HASH(A) returns the SHA-256 hash of the matrix A as a
%   lowercase hexadecimal string. Works for both dense and sparse
%   matrices.
%
%   H = HASH(A, method=method, format=format) computes the hash using the
%   prescribed hash method and output format. 
% 
%   The hash method can be one of {'SHA-1', 'SHA-256', 'MD5'}.
%   The output format can be one of {'hex', 'uint8', 'double'}.

arguments
    A{mustBeNumeric}
    opt.method {mustBeTextScalar} = 'SHA-256'
    opt.format {mustBeTextScalar} = 'hex'
end

%< Convert matrix to uint8 byte stream (preserves sparsity info too)
bytes = getByteStreamFromArray(A);

%< Select algorithm
switch upper(opt.method)
    case 'SHA-256'
        md = java.security.MessageDigest.getInstance('SHA-256');
    case 'SHA-1'
        md = java.security.MessageDigest.getInstance('SHA-1');
    case 'MD5'
        md = java.security.MessageDigest.getInstance('MD5');
    otherwise
        core.except.assert(0, 'UnsupportedMethod', ...
            'method "%s" not supported.', opt.method);
end

md.update(bytes);
raw = typecast(md.digest(), 'uint8');

switch lower(opt.format)
    case 'hex'
        h = lower(dec2hex(raw)')';
        h = reshape(h, 1, []);
    case 'uint8'
        h = raw;
    case 'double'
        h = double(raw);
    otherwise
        core.except.assert(0, 'Unsupportedformat', ...
            'format "%s" not supported.', opt.format);
end
end
