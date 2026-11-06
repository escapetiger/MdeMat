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
%   Supported values are built-in CPU numeric arrays, including sparse and
%   complex arrays. Shape/class/storage are part of identity, as are signed
%   zeros and NaN bit patterns. uint8/double return a row of digest bytes,
%   not a scalar. Corrected hexadecimal identities must be recomputed from
%   arrays; previously permuted hex text is not the canonical digest.

arguments
    A{mustBeNumeric}
    opt.method {mustBeTextScalar} = 'SHA-256'
    opt.format {mustBeTextScalar} = 'hex'
end

if ~ismember(class(A), {'double', 'single', 'int8', 'uint8', 'int16', ...
        'uint16', 'int32', 'uint32', 'int64', 'uint64'})
    error('mdemat:tools:base:hash:UnsupportedType', ...
        'Hashing supports built-in CPU numeric arrays.');
end
if ~ismember(upper(string(opt.method)), ["SHA-256", "SHA-1", "MD5"])
    error('mdemat:tools:base:hash:UnsupportedMethod', ...
        'Hash method "%s" is not supported.', opt.method);
end
if ~ismember(lower(string(opt.format)), ["hex", "uint8", "double"])
    error('mdemat:tools:base:hash:UnsupportedFormat', ...
        'Hash output format "%s" is not supported.', opt.format);
end
if ~usejava('jvm')
    error('mdemat:tools:base:hash:JvmRequired', 'Array hashing requires the MATLAB JVM.');
end

% Canonical CPU-array encoding: class, axes, sparsity, complexity and values.
% This avoids MATLAB's undocumented object serialization. The digest is an
% exact array identity, not a tolerance-based numerical comparison.
header = sprintf('%s|%s|sparse=%d|complex=%d|', class(A), ...
    sprintf('%d,', size(A)), issparse(A), ~isreal(A));
if issparse(A)
    [i, j, values] = find(A);
    payload = [numericBytes(uint64(i)); numericBytes(uint64(j)); ...
        numericBytes(real(values))];
    if ~isreal(A)
        payload = [payload; numericBytes(imag(values))];
    end
else
    payload = numericBytes(real(A));
    if ~isreal(A)
        payload = [payload; numericBytes(imag(A))];
    end
end
bytes = [uint8(header(:)); uint8(0); payload];

%< Select algorithm
switch upper(opt.method)
    case 'SHA-256'
        md = java.security.MessageDigest.getInstance('SHA-256');
    case 'SHA-1'
        md = java.security.MessageDigest.getInstance('SHA-1');
    case 'MD5'
        md = java.security.MessageDigest.getInstance('MD5');
end

md.update(bytes);
raw = reshape(typecast(md.digest(), 'uint8'), 1, []);

switch lower(opt.format)
    case 'hex'
        h = reshape(lower(dec2hex(raw, 2)).', 1, []);
    case 'uint8'
        h = raw;
    case 'double'
        h = double(raw);
end
end

function bytes = numericBytes(values)
[~, ~, endian] = computer;
if endian == 'B'
    values = swapbytes(values);
end
bytes = reshape(typecast(values(:), 'uint8'), [], 1);
end
