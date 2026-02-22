%% EX02_LEGENDRE_BASIS Generate Legendre orthogonal basis functions.
%
%   This script generates Legendre polynomial basis functions up to a
%   specified maximum degree using the improved OrthogonalPolynomialBasis
%   framework. Configure parameters below and run the script.

%% Configuration Parameters
code = 1;
switch code
    case 1
        maxDegree = 10; % Maximum polynomial degree
        variable = sym('x'); % Symbolic variable
        isNormalized = false; % Orthonormal basis flag
        isMonic = true; % Monic polynomials flag
        lowerBound = -1 / 2; % Lower domain bound
        upperBound = 1 / 2; % Upper domain bound
        fileName = '+core/+symbolic/metadata/monic_unit_legendre'; % Save file name
        maxDerivOrder = maxDegree; % Maximum derivative order for compilation
    case 2
        maxDegree = 10; % Maximum polynomial degree
        variable = sym('x'); % Symbolic variable
        isNormalized = true; % Orthonormal basis flag
        isMonic = false; % Monic polynomials flag
        lowerBound = -1; % Lower domain bound
        upperBound = 1; % Upper domain bound
        fileName = '+core/+symbolic/metadata/normal_legendre'; % Save file name
        maxDerivOrder = 0; % Maximum derivative order for compilation
end

%% Generate Legendre Basis Functions

basis = core.symbolic.LegendreBasis(variable, maxDegree, ...
    isNormalized=isNormalized, ...
    isMonic=isMonic, ...
    lower=lowerBound, ...
    upper=upperBound);

basis.compile(maxDerivOrder);

%% Display Results
fprintf('Generated %d basis functions:\n', basis.NCodims);

for i = 1:basis.NCodims
    fprintf('P_%d(x) = %s\n', i-1, char(basis.Expressions(i)));
end
fprintf('\n');

% Display properties
fprintf('Basis Properties:\n');
fprintf('  Max Degree: %d\n', basis.MaxDegree);
fprintf('  IsNormalized: %s\n', string(basis.IsNormalized));
fprintf('  Monic: %s\n', string(basis.IsMonic));
fprintf('  Domain: [%.3f, %.3f]\n', basis.Lower, basis.Upper);
fprintf('  Input Dims: %d\n', basis.NDims);
fprintf('  Output Dims: %d\n', basis.NCodims);
fprintf('\n');

%% Save Results (if specified)
if ~isempty(fileName)
    basis.save(fileName);
    fprintf('Basis saved to: %s\n', [fileName, basis.MetaExt]);
end