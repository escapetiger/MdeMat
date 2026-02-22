%% EX05_FOURIER_BASIS Generate Fourier orthogonal basis functions.
%
%   This script generates trigonometric Fourier basis functions up to a 
%   specified maximum frequency. Configure parameters below and run the 
%   script.

%% Configuration Parameters
maxDegree = 10; % Maximum degree
variable = sym('x'); % Symbolic variable
isNormalized = true; % Orthonormal basis flag
lowerBound = 0; % Lower domain bound
upperBound = 2*pi; % Upper domain bound
fileName = '+core/+symbolic/metadata/normal_fourier'; % Save file name
maxDerivOrder = 0; % Maximum derivative order for compilation

%% Generate Fourier Basis Functions

basis = core.symbolic.FourierBasis(variable, maxDegree, ...
    isNormalized=isNormalized, ...
    lower=lowerBound, ...
    upper=upperBound);
basis.compile(maxDerivOrder);

%% Display Results
fprintf('Generated %d basis functions:\n', basis.NCodims);

for i = 1:basis.NCodims
    fprintf('P_%d(x) = %s\n', i-1, char(basis.Expressions(i)));
end
fprintf('\n');

fprintf('Basis Properties:\n');
fprintf('  Max Degree: %d\n', basis.MaxDegree);
fprintf('  IsNormalized: %s\n', string(basis.IsNormalized));
fprintf('  Domain: [%.3f, %.3f]\n', lowerBound, upperBound);
fprintf('  Input Dims: %d\n', basis.NDims);
fprintf('  Output Dims: %d\n', basis.NCodims);
fprintf('\n');

%% Save Results (if specified)
if ~isempty(fileName)
    basis.save(fileName);
    fprintf('Basis saved to: %s\n', [fileName, basis.MetaExt]);
end