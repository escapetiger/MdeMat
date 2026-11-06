function phi = expseq(Z, n)
%EXPSEQ Evaluate the entire matrix functions phi_0(Z), ..., phi_n(Z).
%   phi{1}=expm(Z) and phi{k+1}=sum_{j>=0} Z^j/(j+k)!. In particular,
%   phi_k(0)=I/k!, including when Z is singular or nilpotent. The identity
%   Z*phi_{k+1}(Z)=phi_k(Z)-I/k! is not used as a linear solve: it does
%   not determine phi_{k+1} on the kernel of Z.
%
%   Small matrices in norm use an adaptive Taylor sum. Otherwise one
%   augmented matrix exponential gives all upper-row blocks without any
%   inverse of Z. Its dimension is (n+1)*size(Z,1); outputs are dense,
%   including for sparse input. Scalar -Inf retains the stiff limit zero
%   for every fixed k. Other nonfinite inputs have no finite matrix value.
%
% See also expm.

if ~(isnumeric(n) && isscalar(n) && isreal(n) && isfinite(n) && n >= 0 && fix(n) == n)
    error('mdemat:ode:InvalidIndex', 'n must be a finite nonnegative integer scalar.');
end
if ~(isnumeric(Z) && ismatrix(Z) && ~isempty(Z))
    error('mdemat:ode:InvalidMatrix', 'Z must be a nonempty numeric matrix.');
end
if size(Z, 1) ~= size(Z, 2)
    error('mdemat:ode:NonSquareMatrix', 'Z must be square.');
end
n = double(n);
if ~isfloat(Z), Z = double(Z); end
Z = full(Z);
phi = cell(1, n+1);
if isscalar(Z) && isreal(Z) && Z == -Inf
    phi(:) = {zeros(size(Z), 'like', Z)};
    return;
end
if any(~isfinite(Z(:)))
    error('mdemat:ode:NonFiniteMatrix', 'Only scalar -Inf is supported as a nonfinite limit.');
end
if n == 0
    phi{1} = expm(Z);
    return;
end
identity = eye(size(Z), 'like', Z);
if norm(Z, 1) <= 0.5
    phi = taylorSequence(Z, identity, n);
    return;
end

% For J=[Z I 0 ...; 0 0 I ...; ...], the first block row of exp(J)
% is [exp(Z),phi_1(Z),...,phi_n(Z)]. Equivalently, phi_k is the exact
% response at t=1 to polynomial forcing t^(k-1)/(k-1)! with zero data.
width = size(Z, 1);
augmented = zeros((n+1)*width, 'like', Z);
augmented(1:width, 1:width) = Z;
for block = 1:n
    rows = (block-1)*width+(1:width);
    augmented(rows, rows+width) = identity;
end
exponential = expm(augmented);
for block = 0:n
    phi{block+1} = exponential(1:width, block*width+(1:width));
end
end

function phi = taylorSequence(Z, identity, n)
% A recurrence for individual Taylor terms avoids powers and cancellation.
phi = cell(1, n+1);
reciprocalFactorial = 1;
for order = 0:n
    if order > 0, reciprocalFactorial = reciprocalFactorial/order; end
    term = reciprocalFactorial*identity;
    value = term;
    for power = 1:128
        term = (term*Z)/(order+power);
        value = value+term;
        if norm(term, 1) <= eps(class(Z))*max(norm(value, 1), realmin(class(Z)))
            break;
        end
    end
    phi{order+1} = value;
end
end
