classdef SeparableFunction < core.linalg.Separable & core.function.Function
    % SEPARABLEFUNCTION Function decomposed as product of factor functions
    % with optional sparse structure.
    %
    %   SeparableFunction represents a mathematical mapping that can be
    %   expressed as a product of separate functions with optional
    %   sparsity:
    %
    %   \f[
    %       f(x_1,x_2,...,x_d) = \sum_{i}\Pi_{j=1}^{d}f_{i_j}(x_j)
    %   \f]
    %
    %   The class supports both dense (full tensor product) and sparse
    %   (reduced set of multi-indices) evaluation patterns. When pattern is
    %   empty ([]), the function uses dense evaluation. When pattern is a
    %   struct with indexer information, it uses sparse evaluation.
    %
    % Examples:
    %   % Dense separable function
    %   factors = {LegendreBasis(3), LegendreBasis(4)};
    %   func = SeparableFunction(factors);
    %
    %   % Sparse separable function with L1 constraint
    %   func = SeparableFunction(factors, 'l1');
    %
    %   % Evaluate function
    %   X = rand(7, 100);  % 7D input (3+4 from factors)
    %   Y = func.evaluate(X);
    %
    % See Also:
    %   core.function.Function, core.linalg.Separable

    properties
        pattern % Sparse pattern specification
    end

    properties (Dependent)
        nFactorCodims % Output dimensions for each factor function
        isSparse % True if function has sparse pattern
    end

    methods
        function obj = SeparableFunction(factors, pattern)
            % SEPARABLEFUNCTION Constructor for separable functions.
            %
            %   obj = SeparableFunction() creates an empty separable
            %   function.
            %
            %   obj = SeparableFunction(factors) creates a separable
            %   function from factor functions using dense evaluation (full
            %   tensor product).
            %
            %   obj = SeparableFunction(factors, pattern) creates a
            %   separable function with specified pattern type for sparse
            %   evaluation.
            %
            % Inputs:
            %   factors - Cell or object array of Function objects (optional)
            %   pattern - Pattern string: {'full' (default), 'l1', 'lx'}
            %
            % Outputs:
            %   obj - SeparableFunction instance

            if nargin < 1, factors = []; end
            if nargin < 2, pattern = 'full'; end

            core.except.assert(ischar(pattern) || isstring(pattern), ...
                'InvalidInput', ...
                'Pattern type must be a string or char array.');

            cls = {'core.function.Function'};
            core.except.assert(core.validate.isAllClass(factors, cls), ...
                'InvalidFactor', ...
                'Factors must be one of {%s}.', strjoin(cls, ', '));

            obj.setFactors(factors);
            obj.setPattern(pattern);
            obj.setNDims();
            obj.setNCodims();
        end

        function obj = setPattern(obj, pattern)
            % SETPATTERN Set the sparse pattern from pattern type.
            %
            %   obj = setPattern(obj, pattern) configures the evaluation
            %   pattern based on the pattern type string. This determines
            %   whether dense or sparse evaluation is used and sets up
            %   appropriate indexing structures.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   pattern - String: 'full', 'l1', or 'lx'
            %
            % Outputs:
            %   obj - The SeparableFunction object

            if isempty(obj.factors)
                obj.pattern = [];
                return;
            end

            pattern = lower(char(pattern));

            supportedTypes = {'full', 'l1', 'lx'};
            core.except.assert( ...
                ismember(pattern, supportedTypes), ...
                'InvalidInput', ...
                'Pattern type must be one of: %s', strjoin(supportedTypes, ', '));

            switch pattern
                case 'full'
                    obj.pattern = [];

                case 'l1'
                    threshold = min(obj.nFactorCodims);

                    indexer = core.linalg.L1MultiIndexer();
                    indexer.setCache(obj.nFactors, obj.nFactors+threshold-1);

                    obj.pattern = struct('indexer', indexer, 'threshold', threshold);

                case 'lx'
                    threshold = min(obj.nFactorCodims);

                    indexer = core.linalg.LxMultiIndexer();
                    indexer.setCache(obj.nFactors, threshold);

                    obj.pattern = struct('indexer', indexer, 'threshold', threshold);
            end
        end

        function Y = evaluate(obj, X)
            % EVALUATE Evaluate the function at specified points.
            %
            %   Y = evaluate(obj, X) evaluates the separable function using
            %   either dense or sparse evaluation based on the pattern
            %   property. Supports both concatenated input format and cell
            %   array input format for factor-wise evaluation.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Coordinates (nDims x ... or cell array)
            %
            % Outputs:
            %   Y - Function values (nCodims x ...)

            core.except.assert(obj.isWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');

            if iscell(X)
                F = obj.factors;
                core.except.assert(length(X) == length(F), ...
                    'InvalidInput', 'Cell inputs mismatch factor count.');

                for i = 1:numel(X)
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end

                    core.except.assert(size(X{i}, 1) == Fi.nDims, ...
                        'InvalidInput', 'Input dimensions mismatch for factor %d.', i);
                end
            else
                n = obj.nDims;
                core.except.assert(size(X, 1) == n, ...
                    'InvalidInput', 'Input dimension mismatch.');
                s = size(X);
                X = reshape(X, n, []);
            end

            Y = obj.evaluateImpl(X);

            m = obj.nCodims;
            core.except.assert(size(Y, 1) == m, ...
                'InvalidOutput', 'Output dimension mismatch.');

            %< reshape output if needed
            if ~iscell(X) && length(s) > 2
                Y = reshape(Y, [m, s(2:end)]);
            end
        end

        function dY = differentiate(obj, X, order)
            % DIFFERENTIATE Compute the derivative of specified order at
            % given points.
            %
            %   dY = differentiate(obj, X, order) computes derivatives of
            %   the separable function using either dense or sparse
            %   evaluation based on the pattern property. Supports both
            %   concatenated input format and cell array format for
            %   factor-wise evaluation.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Coordinates (nDims x ... or cell array)
            %   order - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x ...)

            core.except.assert(obj.isWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');

            n = obj.nDims;
            core.except.assert(isvector(order) && isnumeric(order) ...
                && length(order) == n && all(order >= 0) && sum(order) > 0, ...
                'InvalidInput', 'Bad derivative order.');

            if iscell(X)
                F = obj.factors;
                core.except.assert(length(X) == length(F), ...
                    'InvalidInput', 'Cell inputs mismatch factor count.');

                for i = 1:numel(X)
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end

                    core.except.assert(size(X{i}, 1) == Fi.nDims, ...
                        'InvalidInput', ...
                        'Input dimensions mismatch for factor %d.', i);
                end
            else
                core.except.assert(size(X, 1) == n, ...
                    'InvalidInput', 'Input dimension mismatch.');
                s = size(X);
                X = reshape(X, n, []);
            end

            dY = obj.differentiateImpl(X, order);

            m = obj.nCodims;
            core.except.assert(size(dY, 1) == m, ...
                'InvalidOutput', 'Output dimension mismatch.');

            if ~iscell(X) && length(s) > 2
                dY = reshape(dY, [m, s(2:end)]);
            end
        end

        function n = get.nFactorCodims(obj)
            % GET.NFACTORCODIMS Number of codimensions of each factor.

            F = obj.factors;
            if isempty(F)
                n = [];
                return;
            end

            if iscell(F)
                n = cellfun(@(f) f.nCodims, F);
            else
                n = arrayfun(@(f) f.nCodims, F);
            end
        end

        function TF = get.isSparse(obj)
            % GET.ISSPARSE True if the pattern is non-empty.

            TF = ~isempty(obj.pattern);
        end
    end

    methods (Access = protected)
        function obj = setNDims(obj)
            % SETNDIMS Update input dimensions.
            %
            %   obj = setNDims(obj) computes the total input dimensions as
            %   the sum of input dimensions from all factor functions.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %
            % Outputs:
            %   obj - The SeparableFunction object

            F = obj.factors;
            if isempty(F)
                obj.nDims = 0;
                return;
            end

            if iscell(F)
                obj.nDims = sum(cellfun(@(f) f.nDims, F));
            else
                obj.nDims = sum(arrayfun(@(f) f.nDims, F));
            end
        end

        function obj = setNCodims(obj)
            % SETNCODIMS Update output dimensions.
            %
            %   obj = setNCodims(obj) computes the output dimensions based
            %   on factor functions and pattern specification. For sparse
            %   patterns, nCodims equals the number of active terms in the
            %   indexer cache. For dense patterns, it's the product of all
            %   factor output dimensions.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %
            % Outputs:
            %   obj - The SeparableFunction object

            F = obj.factors;
            if isempty(F)
                obj.nCodims = 0;
                return;
            end

            if obj.isSparse
                obj.nCodims = size(obj.pattern.indexer.cache, 1);
            else
                codims = obj.nFactorCodims;
                obj.nCodims = prod(codims);
            end
        end

        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of function evaluation.
            %
            %   Y = evaluateImpl(obj, X) evaluates the separable function
            %   using either dense or sparse computation based on the
            %   pattern property. This method dispatches to the appropriate
            %   evaluation strategy.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Coordinates (nDims x nPoints or cell array)
            %
            % Outputs:
            %   Y - Function values (nCodims x nPoints)

            if obj.isSparse
                Y = obj.evaluateSparse(X);
            else
                Y = obj.evaluateDense(X);
            end
        end

        function Y = evaluateDense(obj, X)
            % EVALUATEDENSE Dense evaluation using full tensor product.
            %
            %   Y = evaluateDense(obj, X) evaluates the function using the
            %   full tensor product structure. Each factor function is
            %   evaluated and results are combined using Kronecker or
            %   Khatri-Rao products depending on input format.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Coordinates (nDims x nPoints or cell array)
            %
            % Outputs:
            %   Y - Function values (nCodims x nPoints)

            F = obj.factors;
            k = length(F);
            Z = cell(1, k);

            if iscell(X)
                %< evaluate each factor with cell array input
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    mi = Fi.nCodims;
                    Z{i} = Fi.evaluate(X{i});
                    Z{i} = reshape(Z{i}, mi, []);
                end
                Y = core.linalg.kronecker(Z);
            else
                %< evaluate each factor with partitioned input
                j = 1;
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    ni = Fi.nDims;
                    mi = Fi.nCodims;
                    J = j:(j + ni - 1);
                    Z{i} = Fi.evaluate(X(J, :));
                    Z{i} = reshape(Z{i}, mi, []);
                    j = j + ni;
                end
                Y = core.linalg.khatrirao(Z);
            end
        end

        function Y = evaluateSparse(obj, X)
            % EVALUATESPARSE Sparse evaluation using cached multi-indices.
            %
            %   Y = evaluateSparse(obj, X) evaluates the function using a
            %   sparse representation based on cached multi-indices from
            %   the pattern indexer. This reduces computational complexity
            %   for high-dimensional functions.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Coordinates (nDims x nPoints or cell array)
            %
            % Outputs:
            %   Y - Function values (nCodims x nPoints)

            indexer = obj.pattern.indexer;
            M = indexer.cache;

            F = obj.factors;
            k = length(F);

            if iscell(X)
                %< sparse evaluation with cell array input
                Z = cell(1, k);
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    Z{i} = Fi.evaluate(X{i});
                    Z{i} = Z{i}(M(:, i), :);
                end
                Y = prod(cat(3, Z{:}), 3);
            else
                %< sparse evaluation with partitioned input
                Z = cell(1, k);
                j = 1;
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    ni = Fi.nDims;
                    J = j:(j + ni - 1);
                    Z{i} = Fi.evaluate(X(J, :));
                    Z{i} = Z{i}(M(:, i), :);
                    j = j + ni;
                end
                Y = prod(cat(3, Z{:}), 3);
            end
        end

        function dY = differentiateImpl(obj, X, order)
            % DIFFERENTIATEIMPL Implementation of function derivative.
            %
            %   dY = differentiateImpl(obj, X, order) computes derivatives
            %   using either dense or sparse evaluation based on the
            %   pattern property. This method dispatches to the appropriate
            %   evaluation strategy.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Evaluation points (nDims x nPoints or cell array)
            %   order - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            if obj.isSparse
                dY = obj.differentiateSparse(X, order);
            else
                dY = obj.differentiateDense(X, order);
            end
        end

        function dY = differentiateDense(obj, X, order)
            % DIFFERENTIATEDENSE Dense derivative evaluation using full
            % tensor product.
            %
            %   dY = differentiateDense(obj, X, order) computes derivatives using
            %   the full tensor product structure. Each factor function
            %   is evaluated or differentiated as needed, and results are
            %   combined using Kronecker or Khatri-Rao products.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Evaluation points (nDims x nPoints or cell array)
            %   order - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            s = obj.splitOrderVector(order);
            F = obj.factors;
            k = length(F);
            dZ = cell(1, k);

            if iscell(X)
                %< dense derivative evaluation with cell array input
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    mi = Fi.nCodims;
                    if all(s{i} == 0)
                        dZ{i} = Fi.evaluate(X{i});
                    else
                        dZ{i} = Fi.differentiate(X{i}, s{i});
                    end
                    dZ{i} = reshape(dZ{i}, mi, []);
                end
                dY = core.linalg.kronecker(dZ);
            else
                %< dense derivative evaluation with partitioned input
                j = 1;
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    ni = Fi.nDims;
                    mi = Fi.nCodims;
                    J = j:(j + ni - 1);
                    if all(s{i} == 0)
                        dZ{i} = Fi.evaluate(X(J, :));
                    else
                        dZ{i} = Fi.differentiate(X(J, :), s{i});
                    end
                    dZ{i} = reshape(dZ{i}, mi, []);
                    j = j + ni;
                end
                dY = core.linalg.khatrirao(dZ);
            end
        end

        function dY = differentiateSparse(obj, X, order)
            % DIFFERENTIATESPARSE Sparse derivative evaluation using cached
            % multi-indices.
            %
            %   dY = differentiateSparse(obj, X, order) computes
            %   derivatives using a sparse representation based on cached
            %   multi-indices from the pattern indexer. This reduces
            %   computational complexity for high-dimensional functions.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   X - Evaluation points (nDims x nPoints or cell array)
            %   order - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            s = obj.splitOrderVector(order);
            indexer = obj.pattern.indexer;
            M = indexer.cache;

            F = obj.factors;
            k = length(F);

            if iscell(X)
                %< sparse derivative evaluation with cell array input
                dZ = cell(1, k);
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    if all(s{i} == 0)
                        dZ{i} = Fi.evaluate(X{i});
                    else
                        dZ{i} = Fi.differentiate(X{i}, s{i});
                    end
                    dZ{i} = dZ{i}(M(:, i), :);
                end
                dY = prod(cat(3, dZ{:}), 3);
            else
                %< sparse derivative evaluation with partitioned input
                dZ = cell(1, k);
                j = 1;
                for i = 1:k
                    if iscell(F)
                        Fi = F{i};
                    else
                        Fi = F(i);
                    end
                    ni = Fi.nDims;
                    J = j:(j + ni - 1);
                    if all(s{i} == 0)
                        dZ{i} = Fi.evaluate(X(J, :));
                    else
                        dZ{i} = Fi.differentiate(X(J, :), s{i});
                    end
                    dZ{i} = dZ{i}(M(:, i), :);
                    j = j + ni;
                end
                dY = prod(cat(3, dZ{:}), 3);
            end
        end

        function s = splitOrderVector(obj, order)
            % SPLITORDERVECTOR Split derivative order vector for each
            % factor.
            %
            %   s = splitOrderVector(obj, order) decomposes a global
            %   derivative order vector into separate order vectors for
            %   each factor function, enabling factor-wise derivative
            %   computation.
            %
            % Inputs:
            %   obj - The SeparableFunction object
            %   order - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   s - Cell array of order vectors for each factor

            F = obj.factors;
            k = length(F);
            s = cell(1, k);
            j = 1;
            for i = 1:k
                if iscell(F)
                    Fi = F{i};
                else
                    Fi = F(i);
                end
                ni = Fi.nDims;
                s{i} = order(j:(j + ni - 1));
                j = j + ni;
            end
        end
    end
end