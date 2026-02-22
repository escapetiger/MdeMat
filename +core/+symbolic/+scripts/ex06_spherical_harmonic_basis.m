%% EX06_SPHERICAL_HARMONIC_BASIS Generate spherical harmonic orthogonal basis functions.
%
%   This script generates spherical harmonic basis functions up to a
%   specified maximum degree. Configure parameters below and run the
%   script.

%% Configuration Parameters
maxDegree = 5; % Maximum degree l
variables = [sym('mu'), sym('theta')]; % Symbolic variables (mu, theta)
isNormalized = true; % Orthonormal basis flag
fileName = '+core/+symbolic/metadata/normal_spherical_harmonic'; % Save file name
maxDerivOrder = 0; % Maximum derivative order for handles

%% Generate Spherical Harmonic Basis Functions
fprintf('=== Spherical Harmonic Basis Generation ===\n');
fprintf('Generating spherical harmonic basis up to degree %d...\n', ...
    maxDegree);

basis = core.symbolic.SphericalHarmonicBasis(variables, maxDegree, ...
    isNormalized=isNormalized);
basis.compile(maxDerivOrder);

%% Display Results
fprintf('Generated %d basis functions:\n', basis.NCodims);

for i = 1:basis.NCodims
    fprintf('P_%d(mu,theta) = %s\n', i-1, char(basis.Expressions(i)));
end
fprintf('\n');

fprintf('Basis Properties:\n');
fprintf('  Max Degree: %d\n', basis.MaxDegree);
fprintf('  IsNormalized: %s\n', string(basis.IsNormalized));
fprintf('  Input Dims: %d\n', basis.NDims);
fprintf('  Output Dims: %d\n', basis.NCodims);
fprintf('\n');

%% Save Results (if specified)
if ~isempty(fileName)
    basis.save(fileName);
    fprintf('Basis saved to: %s\n', [fileName, basis.MetaExt]);
end