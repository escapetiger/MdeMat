function ex01_lagrange_basis(options)
%EX01_LAGRANGE_BASIS Regenerate Gauss interpolation metadata explicitly.
%   Gauss-Legendre six-point roots use the displayed 34-digit approximations
%   to the roots of P6; the resulting runtime basis is evaluated in double.
arguments
    options.overwrite (1, 1) logical = false
end
%
%   This script generates Lagrange basis functions for specified
%   interpolation nodes. Configure parameters below and run the script.

%% Configuration Parameters
metadataDirectory = fullfile(fileparts(mfilename('fullpath')), 'metadata');

assert(options.overwrite, 'MdeMat:BasisGeneration:OverwriteRequired', ...
    'Regeneration replaces development metadata; set overwrite=true explicitly.');
for code = 1:2
switch code
    case 1
        nodes = { ...
            [sym(0)]; ... % n=1
            [-sym(1) / (2 * sym(sqrt(3))), ...
            sym(1) / (2 * sym(sqrt(3)))]; ... % n=2
            [sym(0), ...
            -sym(sqrt(sym(3)/5)) / 2, ...
            sym(sqrt(sym(3)/5)) / 2]; ... % n=3
            [-sqrt((15+2*sqrt(sym(30)))/35)/2, ...
            sqrt((15+2*sqrt(sym(30)))/35)/2, ...
            -sqrt((15-2*sqrt(sym(30)))/35)/2, ...
            sqrt((15-2*sqrt(sym(30)))/35)/2]; ... % n=4
            [sym(0), ...
            -sym(1) / 6 * sym(sqrt(5-2*sqrt(sym(10)/7))), ...
            sym(1) / 6 * sym(sqrt(5-2*sqrt(sym(10)/7))), ...
            -sym(1) / 6 * sym(sqrt(5+2*sqrt(sym(10)/7))), ...
            sym(1) / 6 * sym(sqrt(5+2*sqrt(sym(10)/7)))]; ... % n=5
            [-sym('0.9324695142031520278123015544939946'), ...
            sym('0.9324695142031520278123015544939946'), ...
            -sym('0.6612093864662645136613995950199053'), ...
            sym('0.6612093864662645136613995950199053'), ...
            -sym('0.2386191860831969086305017216807119'), ...
            sym('0.2386191860831969086305017216807119')]/2; ... % n=6
            };
        fileName = fullfile(metadataDirectory, 'unit_gauss_legendre_lagrange');
    case 2
        nodes = { ...
            [sym(0)]; ... % n=1
            [-sym(1) / 2, sym(1) / 2]; ... % n=2
            [-sym(1) / 2, sym(0), sym(1) / 2]; ... % n=3
            [-sym(1) / 2, -sym(sqrt(sym(1)/5)) / 2, sym(sqrt(sym(1)/5)) / 2, sym(1) / 2]; ... % n=4
            [-sym(1) / 2, -sym(sqrt(sym(3)/7)) / 2, sym(0), sym(sqrt(sym(3)/7)) / 2, sym(1) / 2]; ... % n=5
            [-sym(1)/2, -sqrt((7+2*sqrt(sym(7)))/21)/2, ...
            -sqrt((7-2*sqrt(sym(7)))/21)/2, sqrt((7-2*sqrt(sym(7)))/21)/2, ...
            sqrt((7+2*sqrt(sym(7)))/21)/2, sym(1)/2]; ... % n=6
            };
        fileName = fullfile(metadataDirectory, 'unit_gauss_lobatto_lagrange');
end

variable = sym('x'); % Symbolic variable
maxDerivOrder = length(nodes) - 1; % Maximum derivative order for handles

%% Generate Lagrange Basis Functions

for k = 1:length(nodes)
    basis = mdemat.common.symbolic.LagrangeBasis(variable, nodes{k});

    basis.compile(maxDerivOrder);

    fprintf('Generated %d basis functions:\n', basis.NCodims);

    for i = 1:basis.NCodims
        fprintf('P_%d(x) = %s\n', i-1, char(basis.Expressions(i)));
    end
    fprintf('\n');

    fprintf('Basis Properties:\n');
    fprintf('  Nodes: %s\n', string(basis.Nodes));
    fprintf('  Input Dims: %d\n', basis.NDims);
    fprintf('  Output Dims: %d\n', basis.NCodims);
    fprintf('\n');

    if ~isempty(fileName)
        basis.save(sprintf('%s_%d', fileName, k));
        fprintf('Basis saved to: %s\n', [sprintf('%s_%d', fileName, k), basis.MetaExt]);
    end
end
end
end
