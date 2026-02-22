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
    % See Also:
    %   core.function.Function, core.linalg.Separable
    
    properties
        Pattern % Sparse pattern specification
    end
    
    properties (Dependent)
        NFactorCodims % Output dimensions for each factor function
        IsSparse % True if function has sparse pattern
    end
    
    methods
        function obj = SeparableFunction(options)
            % SEPARABLEFUNCTION Constructor for separable functions.
            %
            %   obj = SeparableFunction() creates an empty separable
            %   function.
            %
            %   obj = SeparableFunction('factors', factors) creates a
            %   separable function from factor functions using dense
            %   evaluation.
            %
            %   obj = SeparableFunction('factors', factors, 'pattern', pattern)
            %   creates a separable function with specified pattern type for
            %   sparse evaluation.
            
            arguments
                options.factors{core.except.mustBeCellOrObject} = {}
                options.pattern{mustBeMember(options.pattern, {'', 'P', 'Q'})} = ''
            end
            
            cls = {'core.function.Function'};
            isValidFactors = core.except.isAllClass(options.factors, cls);
            core.except.assert(isValidFactors, 'InvalidFactor', ...
                'Factors must be one of {%s}.', strjoin(cls, ', '));
            
            obj.setFactors(options.factors);
            obj.setPattern(options.pattern);
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
            
            arguments
                obj core.function.SeparableFunction
                pattern{mustBeMember(pattern, {'', 'P', 'Q'})} = ''
            end
            
            if isempty(obj.Factors)
                obj.Pattern = [];
                return;
            end
            
            switch pattern
                case 'P'
                    threshold = min(obj.NFactorCodims);
                    
                    indexer = core.linalg.L1MultiIndexer();
                    indexer.setCache(obj.NFactors, obj.NFactors+threshold-1);
                    
                    obj.Pattern = struct('indexer', indexer, ...
                        'threshold', threshold);
                    
                case 'Q'
                    threshold = min(obj.NFactorCodims);
                    
                    indexer = core.linalg.LxMultiIndexer();
                    indexer.setCache(obj.NFactors, threshold);
                    
                    obj.Pattern = struct('indexer', indexer, ...
                        'threshold', threshold);
                    
                otherwise
                    obj.Pattern = [];
            end
        end
        
        function Y = eval(obj, X)
            % EVAL Evaluate the function at specified points.
            %
            %   Y = eval(obj, X) evaluates the separable function using
            %   either dense or sparse evaluation based on the pattern
            %   property. Supports both concatenated input format and cell
            %   array input format for factor-wise evaluation.
            
            arguments
                obj core.function.SeparableFunction
                X{mustBeA(X, {'numeric', 'cell'})}
            end
            
            core.except.assert(obj.IsWellDefined, 'NotWellDefined', ...
                'Function is not well-defined.');
            
            if iscell(X)
                isValidCellInput = length(X) == obj.NFactors;
                core.except.assert(isValidCellInput, 'InvalidInput', ...
                    'Cell inputs mismatch factor count.');
                
                for i = 1:numel(X)
                    F = obj.getFactor(i);
                    isValidInput = size(X{i}, 1) == F.NDims;
                    core.except.assert(isValidInput, 'InvalidInput', ...
                        'Input dimensions mismatch for factor %d.', i);
                end
            else
                nd = obj.NDims;
                isValidNumericInput = size(X, 1) == nd;
                core.except.assert(isValidNumericInput, 'InvalidInput', ...
                    'Input dimension mismatch.');
                nx = size(X);
                X = reshape(X, nd, []);
            end
            
            Y = obj.evalImpl(X);
            
            nc = obj.NCodims;
            isValidOutput = size(Y, 1) == nc;
            core.except.assert(isValidOutput, 'InvalidOutput', ...
                'Output dimension mismatch.');
            
            if ~iscell(X) && length(nx) > 2
                Y = reshape(Y, [nc, nx(2:end)]);
            end
        end
        
        function dY = diff(obj, X, order)
            % DIFF Compute the derivative of specified order at
            % given points.
            %
            %   dY = diff(obj, X, order) computes derivatives of
            %   the separable function using either dense or sparse
            %   evaluation based on the pattern property. Supports both
            %   concatenated input format and cell array format for
            %   factor-wise evaluation.
            
            arguments
                obj core.function.SeparableFunction
                X{mustBeA(X, {'numeric', 'cell'})}
                order{mustBeVector, mustBeNonnegative}
            end
            
            core.except.assert(obj.IsWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');
            
            nd = obj.NDims;
            isValidOrder = length(order) == nd && sum(order) > 0;
            core.except.assert(isValidOrder, 'InvalidInput', ...
                'Bad derivative order.');
            
            if iscell(X)
                isValidCellInput = length(X) == obj.NFactors;
                core.except.assert(isValidCellInput, 'InvalidInput', ...
                    'Cell inputs mismatch factor count.');
                
                for i = 1:numel(X)
                    F = obj.getFactor(i);
                    isValidInput = size(X{i}, 1) == F.NDims;
                    core.except.assert(isValidInput, 'InvalidInput', ...
                        'Input dimensions mismatch for factor %d.', i);
                end
            else
                nd = obj.NDims;
                isValidNumericInput = size(X, 1) == nd;
                core.except.assert(isValidNumericInput, 'InvalidInput', ...
                    'Input dimension mismatch.');
                nx = size(X);
                X = reshape(X, nd, []);
            end
            
            dY = obj.diffImpl(X, order);
            
            nc = obj.NCodims;
            isValidOutput = size(dY, 1) == nc;
            core.except.assert(isValidOutput, 'InvalidOutput', ...
                'Output dimension mismatch.');
            
            if ~iscell(X) && length(nx) > 2
                dY = reshape(dY, [nc, nx(2:end)]);
            end
        end
        
        function n = get.NFactorCodims(obj)
            % GET.NFACTORCODIMS Number of codimensions of each factor.
            
            F = obj.Factors;
            if isempty(F)
                n = [];
                return;
            end
            
            if iscell(F)
                n = cellfun(@(f) f.NCodims, F);
            else
                n = arrayfun(@(f) f.NCodims, F);
            end
        end
        
        function TF = get.IsSparse(obj)
            % GET.ISSPARSE True if the pattern is non-empty.
            
            TF = ~isempty(obj.Pattern);
        end
    end
    
    methods (Access = protected)
        function obj = setNDims(obj)
            % SETNDIMS Update input dimensions.
            %
            %   obj = setNDims(obj) computes the total input dimensions as
            %   the sum of input dimensions from all factor functions.
            
            F = obj.Factors;
            if isempty(F)
                obj.NDims = 0;
                return;
            end
            
            if iscell(F)
                obj.NDims = sum(cellfun(@(f) f.NDims, F));
            else
                obj.NDims = sum(arrayfun(@(f) f.NDims, F));
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
            
            F = obj.Factors;
            if isempty(F)
                obj.NCodims = 0;
                return;
            end
            
            if obj.IsSparse
                obj.NCodims = size(obj.Pattern.indexer.Cache, 1);
            else
                obj.NCodims = prod(obj.NFactorCodims);
            end
        end
        
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.
            %
            %   Y = evalImpl(obj, X) evaluates the separable function
            %   using either dense or sparse computation based on the
            %   pattern property. This method dispatches to the appropriate
            %   evaluation strategy.
            
            if obj.IsSparse
                Y = obj.evalSparse(X);
            else
                Y = obj.evalDense(X);
            end
        end
        
        function Y = evalDense(obj, X)
            % EVALDENSE Dense evaluation using full tensor product.
            %
            %   Y = evalDense(obj, X) evaluates the function using the
            %   full tensor product structure. Each factor function is
            %   evaluated and results are combined using Kronecker or
            %   Khatri-Rao products depending on input format.
            
            nf = obj.NFactors;
            Z = cell(1, nf);
            
            if iscell(X)
                %< evaluate each factor with cell array input
                for i = 1:nf
                    F = obj.getFactor(i);
                    nc = F.NCodims;
                    Z{i} = F.eval(X{i});
                    Z{i} = reshape(Z{i}, nc, []);
                end
                Y = core.linalg.kronecker(Z);
            else
                %< evaluate each factor with partitioned input
                j = 1;
                for i = 1:nf
                    F = obj.getFactor(i);
                    nd = F.NDims;
                    nc = F.NCodims;
                    J = j:(j + nd - 1);
                    Z{i} = F.eval(X(J, :));
                    Z{i} = reshape(Z{i}, nc, []);
                    j = j + nd;
                end
                Y = core.linalg.khatrirao(Z);
            end
        end
        
        function Y = evalSparse(obj, X)
            % EVALSPARSE Sparse evaluation using cached multi-indices.
            %
            %   Y = evalSparse(obj, X) evaluates the function using a
            %   sparse representation based on cached multi-indices from
            %   the pattern indexer. This reduces computational complexity
            %   for high-dimensional functions.
            
            indexer = obj.Pattern.indexer;
            M = indexer.Cache;
            nf = obj.NFactors;
            
            Z = cell(1, nf);
            
            if iscell(X)
                %< sparse evaluation with cell array input
                for i = 1:nf
                    F = obj.getFactor(i);
                    Z{i} = F.eval(X{i});
                    Z{i} = Z{i}(M(:, i), :);
                end
                Y = prod(cat(3, Z{:}), 3);
            else
                %< sparse evaluation with partitioned input
                j = 1;
                for i = 1:nf
                    F = obj.getFactor(i);
                    nd = F.NDims;
                    J = j:(j + nd - 1);
                    Z{i} = F.eval(X(J, :));
                    Z{i} = Z{i}(M(:, i), :);
                    j = j + nd;
                end
                Y = prod(cat(3, Z{:}), 3);
            end
        end
        
        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of function derivative.
            %
            %   dY = diffImpl(obj, X, order) computes derivatives
            %   using either dense or sparse evaluation based on the
            %   pattern property. This method dispatches to the appropriate
            %   evaluation strategy.
            
            if obj.IsSparse
                dY = obj.diffSparse(X, order);
            else
                dY = obj.diffDense(X, order);
            end
        end
        
        function dY = diffDense(obj, X, order)
            % DIFFDENSE Dense derivative evaluation using full
            % tensor product.
            %
            %   dY = diffDense(obj, X, order) computes derivatives using
            %   the full tensor product structure. Each factor function
            %   is evaluated or differentiated as needed, and results are
            %   combined using Kronecker or Khatri-Rao products.
            
            m = obj.splitOrderVector(order);
            nf = obj.NFactors;
            dZ = cell(1, nf);
            
            if iscell(X)
                %< dense derivative evaluation with cell array input
                for i = 1:k
                    F = obj.getFactor(i);
                    nc = F.NCodims;
                    if all(m{i} == 0)
                        dZ{i} = F.eval(X{i});
                    else
                        dZ{i} = F.diff(X{i}, m{i});
                    end
                    dZ{i} = reshape(dZ{i}, nc, []);
                end
                dY = core.linalg.kronecker(dZ);
            else
                %< dense derivative evaluation with partitioned input
                j = 1;
                for i = 1:nf
                    F = obj.getFactor(i);
                    nd = F.NDims;
                    nc = F.NCodims;
                    J = j:(j + nd - 1);
                    if all(m{i} == 0)
                        dZ{i} = F.eval(X(J, :));
                    else
                        dZ{i} = F.diff(X(J, :), m{i});
                    end
                    dZ{i} = reshape(dZ{i}, nc, []);
                    j = j + nd;
                end
                dY = core.linalg.khatrirao(dZ);
            end
        end
        
        function dY = diffSparse(obj, X, order)
            % DIFFSPARSE Sparse derivative evaluation using cached
            % multi-indices.
            %
            %   dY = diffSparse(obj, X, order) computes derivatives using a
            %   sparse representation based on cached multi-indices from
            %   the pattern indexer. This reduces computational complexity
            %   for high-dimensional functions.
            
            m = obj.splitOrderVector(order);
            nf = obj.NFactors;
            dZ = cell(1, nf);
            indexer = obj.Pattern.indexer;
            M = indexer.Cache;
            
            if iscell(X)
                %< sparse derivative evaluation with cell array input
                for i = 1:k
                    F = obj.getFactor(i);
                    if all(m{i} == 0)
                        dZ{i} = F.eval(X{i});
                    else
                        dZ{i} = F.diff(X{i}, m{i});
                    end
                    dZ{i} = dZ{i}(M(:, i), :);
                end
                dY = prod(cat(3, dZ{:}), 3);
            else
                %< sparse derivative evaluation with partitioned input
                j = 1;
                for i = 1:nf
                    F = obj.getFactor(i);
                    nd = F.NDims;
                    J = j:(j + nd - 1);
                    if all(m{i} == 0)
                        dZ{i} = F.eval(X(J, :));
                    else
                        dZ{i} = F.diff(X(J, :), m{i});
                    end
                    dZ{i} = dZ{i}(M(:, i), :);
                    j = j + nd;
                end
                dY = prod(cat(3, dZ{:}), 3);
            end
        end
        
        function m = splitOrderVector(obj, order)
            % SPLITORDERVECTOR Split derivative order vector for each
            % factor.
            %
            %   m = splitOrderVector(obj, order) decomposes a global
            %   derivative order vector into separate order vectors for
            %   each factor function, enabling factor-wise derivative
            %   computation.
            
            nf = obj.NFactors;
            m = cell(1, nf);
            j = 1;
            for i = 1:nf
                F = obj.getFactor(i);
                nd = F.NDims;
                m{i} = order(j:(j + nd - 1));
                j = j + nd;
            end
        end
    end
end