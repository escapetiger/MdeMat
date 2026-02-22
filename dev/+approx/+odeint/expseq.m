function phi = expseq(Z, n)
    % EXPSEQ Generate exponential sequence phi_k(Z).
    %
    %   phi = expseq(Z, n) computes the first (n+1) terms of the
    %   exponential sequence defined by:
    %
    %   \f[
    %     \phi_0(Z) = \exp(Z)
    %   \f]
    %   \f[
    %     \phi_{k+1}(Z) = \frac{\phi_k(Z) - \phi_k(0)}{Z}
    %   \f]
    %   \f[
    %     \phi_k(0) = \frac{1}{k!}
    %   \f]
    %
    %   The sequence is fundamental in exponential integrators and provides
    %   exact integration of polynomial source terms multiplied by
    %   exponential decay. These functions are also known as the matrix
    %   \f$\phi\f$-functions and appear in the exact solution of linear
    %   ODEs with polynomial forcing terms.
    %
    %   For small ||Z||, series expansions are used to avoid numerical
    %   issues with division by Z. The function automatically handles
    %   both sparse and dense matrices appropriately. The implementation
    %   switches between series expansion and recurrence relation based
    %   on the Frobenius norm of Z.
    %
    % See also:
    %   expm, factorial, approx.odeint.ExponentialIntegrator

    %< Validate inputs
    core.except.assert(n >= 0, ...
        'InvalidIndex', 'Index n must be non-negative.');
    core.except.assert(ismatrix(Z) && isnumeric(Z), ...
        'InvalidMatrix', 'Z must be a numeric matrix.');
    core.except.assert(size(Z, 1) == size(Z, 2), ...
        'NonSquareMatrix', 'Z must be a square matrix.');

    %< Initialize output
    phi = cell(1, n + 1);
    
    %< Determine matrix type for identity matrix
    if issparse(Z)
        I = speye(size(Z));
    else
        I = eye(size(Z));
    end
    
    %< Compute phi_0(Z) = exp(Z)
    if isinf(Z) && all(Z < 0)
        phi{1} = zeros(size(Z));
    else
        phi{1} = expm(Z);
    end
    
    %< Return if only phi_0 is requested
    if n == 0
        return;
    end
    
    %< Check if we should use series expansion for numerical stability
    useSeriesExpansion = norm(Z, 'fro') < 1e-2;
    
    if useSeriesExpansion
        %< Use series expansion for all terms when Z is small
        phiTerms = computeSeriesExpansion(Z, I, n);
        for k = 1:n
            phi{k + 1} = phiTerms{k};
        end
    else
        %< Use recurrence relation for larger Z
        phiTerms = computeRecurrenceRelation(Z, I, phi{1}, n);
        for k = 1:n
            phi{k + 1} = phiTerms{k};
        end
    end
end

function phiTerms = computeSeriesExpansion(Z, I, n)
    % COMPUTESERIESEXPANSION Compute phi_k using series expansion.
    %
    %   phiTerms = computeSeriesExpansion(Z, I, n) computes phi_k(Z) using
    %   the Taylor series expansion:
    %
    %   \f[
    %     \phi_k(Z) = \sum_{j=0}^{\infty} \frac{Z^j}{(k+j)!}
    %   \f]
    %
    %   This series expansion is used for numerical stability when ||Z|| is
    %   small, avoiding the division by Z in the recurrence relation.

    phiTerms = cell(n, 1);
    
    %< Pre-compute powers of Z for efficiency
    maxTerms = 20; %< Sufficient for convergence when ||Z|| < 1e-2
    ZPowers = cell(maxTerms, 1);
    ZPowers{1} = I; %< Z^0
    
    for j = 2:maxTerms
        ZPowers{j} = ZPowers{j-1} * Z;
    end
    
    %< Compute each phi_k using series expansion
    for k = 1:n
        phiK = zeros(size(Z));
        
        %< Sum series: phi_k(Z) = sum_{j=0}^{inf} Z^j / (k+j)!
        for j = 1:maxTerms
            coefficient = 1 / factorial(k + j - 1);
            term = coefficient * ZPowers{j};
            phiK = phiK + term;
            
            %< Check convergence
            if j > 5 && norm(term, 'fro') < 1e-15 * norm(phiK, 'fro')
                break;
            end
        end
        
        phiTerms{k} = phiK;
    end
end

function phiTerms = computeRecurrenceRelation(Z, I, phi0, n)
    % COMPUTERECURRENCERELATION Compute phi_k using recurrence relation.
    %
    %   phiTerms = computeRecurrenceRelation(Z, I, phi0, n) computes
    %   phi_{k+1}(Z) = (phi_k(Z) - phi_k(0)) / Z for k = 0, ..., n-1
    %   using the recurrence relation for numerical stability when ||Z|| is
    %   not small.

    phiTerms = cell(n, 1);
    
    %< Compute phi_k(0) = 1/k! for all needed values
    phiAtZero = zeros(n + 1, 1);
    for k = 0:n
        phiAtZero(k + 1) = 1 / factorial(k);
    end
    
    %< Initialize with phi_0
    currentPhi = phi0;
    
    %< Apply recurrence relation
    for k = 1:n
        %< phi_{k+1}(Z) = (phi_k(Z) - phi_k(0)) / Z
        numerator = currentPhi - phiAtZero(k) * I;
        
        %< Solve Z * phi_{k+1} = numerator for phi_{k+1}
        %< This avoids explicit matrix inversion
        phiTerms{k} = Z \ numerator;
        
        %< Update for next iteration
        currentPhi = phiTerms{k};
    end
end