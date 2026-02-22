%% EX01_LAGRANGE_BASIS Generate Lagrange interpolation basis functions.
%
%   This script generates Lagrange basis functions for specified
%   interpolation nodes. Configure parameters below and run the script.

%% Configuration Parameters

code = 2;
switch code
    case 1
        nodes = { ...
            [sym(0)]; ... % n=1
            [-sym(1) / (2 * sym(sqrt(3))), ...
            sym(1) / (2 * sym(sqrt(3)))]; ... % n=2
            [sym(0), ...
            -sym(sqrt(sym(3)/5)) / 2, ...
            sym(sqrt(sym(3)/5)) / 2]; ... % n=3
            [-sym(1) / 70 * sym(sqrt(30+2*sqrt(sym(30)))), ...
            sym(1) / 70 * sym(sqrt(30+2*sqrt(sym(30)))), ...
            -sym(1) / 70 * sym(sqrt(30-2*sqrt(sym(30)))), ...
            sym(1) / 70 * sym(sqrt(30-2*sqrt(sym(30))))]; ... % n=4
            [sym(0), ...
            -sym(1) / 6 * sym(sqrt(5-2*sqrt(sym(10)/7))), ...
            sym(1) / 6 * sym(sqrt(5-2*sqrt(sym(10)/7))), ...
            -sym(1) / 6 * sym(sqrt(5+2*sqrt(sym(10)/7))), ...
            sym(1) / 6 * sym(sqrt(5+2*sqrt(sym(10)/7)))]; ... % n=5
            [-sym(1) / 42 * sym(sqrt(525+70*sqrt(sym(30)))), ...
            sym(1) / 42 * sym(sqrt(525+70*sqrt(sym(30)))), ...
            -sym(1) / 42 * sym(sqrt(525-70*sqrt(sym(30)))), ...
            sym(1) / 42 * sym(sqrt(525-70*sqrt(sym(30)))), ...
            -sym(1) / 14 * sym(sqrt(5)), ...
            sym(1) / 14 * sym(sqrt(5))]; ... % n=6
            };
        fileName = '+core/+symbolic/metadata/unit_gauss_legendre_lagrange'; % Save file name
    case 2
        nodes = { ...
            [sym(0)]; ... % n=1
            [-sym(1) / 2, sym(1) / 2]; ... % n=2
            [-sym(1) / 2, sym(0), sym(1) / 2]; ... % n=3
            [-sym(1) / 2, -sym(sqrt(sym(1)/5)) / 2, sym(sqrt(sym(1)/5)) / 2, sym(1) / 2]; ... % n=4
            [-sym(1) / 2, -sym(sqrt(sym(3)/7)) / 2, sym(0), sym(sqrt(sym(3)/7)) / 2, sym(1) / 2]; ... % n=5
            [-sym(1) / 2, -sym(sqrt(sym(1)/3)) / 2, -sym(sqrt(sym(3)/7)) / 2, ...
            sym(sqrt(sym(3)/7)) / 2, sym(sqrt(sym(1)/3)) / 2, sym(1) / 2]; ... % n=6
            };
        fileName = '+core/+symbolic/metadata/unit_gauss_lobatto_lagrange'; % Save file name
end

variable = sym('x'); % Symbolic variable
maxDerivOrder = length(nodes) - 1; % Maximum derivative order for handles

%% Generate Lagrange Basis Functions

for k = 1:length(nodes)
    basis = core.symbolic.LagrangeBasis(variable, nodes{k});

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