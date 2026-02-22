%% EX04_HERMITE_BASIS Generate probabilist's Hermite orthogonal basis functions.
%
%   This script generates probabilist's Hermite polynomial basis functions 
%   up to a specified maximum degree. Configure parameters below and run 
%   the script.

%% Configuration Parameters
maxDegree = 10; % Maximum polynomial degree
variable = sym('x'); % Symbolic variable
isNormalized = true; % Orthonormal basis flag
isMonic = false; % Monic polynomials flag
fileName = '+core/+symbolic/metadata/normal_hermite'; % Save file name
maxDerivOrder = maxDegree; % Maximum derivative order for compilation

%% Generate Hermite Basis Functions

basis = core.symbolic.HermiteBasis(variable, maxDegree, ...
    isNormalized=isNormalized, ...
    isMonic=isMonic);

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
fprintf('  Domain: [%.3f, %.3f]\n', -inf, inf);
fprintf('  Input Dims: %d\n', basis.NDims);
fprintf('  Output Dims: %d\n', basis.NCodims);
fprintf('\n');

%% Save Results (if specified)
if ~isempty(fileName)
    basis.save(fileName);
    fprintf('Basis saved to: %s\n', [fileName, basis.MetaExt]);
end