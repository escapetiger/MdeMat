classdef SumSpace < core.linalg.Separable
    % SUMSPACE Direct sum of spectral function spaces.
    %
    %   SumSpace represents the direct sum of spectral spaces of which each
    %   entry has the form:
    %
    %   \f[
    %       f(x) = \sum_{i=1}^k \sum_{j=1}^{n_i} w_{i,j} c_{i,j}
    %       \phi_{i,j}(x)
    %   \f]
    %
    %   where \f$\phi_{i,j}(x)\f$ are basis functions from space \f$i\f$,
    %   \f$c_{i,j}\f$ are coefficients, and \f$w_{i,j}\f$ are per-DOF
    %   weights.
    %
    %   The direct sum space \f$W = V_1 \oplus V_2 \oplus \cdots \oplus
    %   V_k\f$ has dimension \f$\dim(W) = \sum_i \dim(V_i)\f$.
    %
    % Examples:
    %   % Create direct sum from multiple spectral spaces
    %   space1 = approx.space.SpectralSpace(element1);
    %   space2 = approx.space.SpectralSpace(element2);
    %   weights1 = ones(space1.nDofs, 1);      % Per-DOF weights
    %   weights2 = 2 * ones(space2.nDofs, 1);  % Scaled weights
    %   sumSpace = approx.space.SumSpace({space1, space2}, {weights1, weights2});
    %
    %   % Evaluate with separate coefficient matrices
    %   points = [0, 0.5, 1.0; 0, 0, 0];  % 2D points
    %   coeffs1 = [1; 0.5];               % Coefficients for space1
    %   coeffs2 = [-0.2; 0.3; 0.1];       % Coefficients for space2
    %   values = sumSpace.evaluate(points, coeffs1, coeffs2);
    %
    % See also:
    %   approx.space.SpectralSpace

    properties
        weights % Cell array of weight vectors, one per space
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom (integer)
        nDims % Number of spatial dimensions (integer)
    end

    methods
        function obj = SumSpace(factors, weights)
            % SUMSPACE Constructor for SumSpace.
            %
            %   obj = SumSpace(factors, weights) creates a direct sum
            %   of the specified spectral spaces with the given per-DOF
            %   weights. All spaces must have the same spatial dimension.
            %
            % Inputs:
            %   factors - Cell array of SpectralSpace objects
            %   weights - Cell array of weight vectors (optional)
            %
            % Outputs:
            %   obj - Constructed SumSpace object

            cls = {'approx.space.SpectralSpace'};
            core.except.assert( ...
                core.validate.isAllClass(factors, cls), ...
                'InvalidFactor', ...
                'Factors must be one of {%s}.', strjoin(cls, ', '));

            obj.setFactors(factors);
            obj.setWeights(weights);
        end

        function obj = setWeights(obj, weights)
            % SETWEIGHTS Set the weights.
            %
            %   obj = setWeights(obj, weights) sets weights @a weights for
            %   the direct sum object.
            %
            % Inputs:
            %   obj - The SumSpace object
            %   weights - Factor weights (cell array) 
            %
            % Outputs:
            %   obj - The SumSpace object

            if isempty(weights)
                obj.weights = cell(1, obj.nFactors);
                return; 
            end

            core.except.assert(length(weights) == obj.nFactors, ...
                'InvalidInput', ...
                'Number of weight vectors must match number of spaces.');

            tf = arrayfun(@(i) length(weights{i}) == obj.getFactor(i).nDofs, 1:obj.nFactors);
            core.except.assert(all(tf), 'InvalidInput', ...
                'Weight vector must match the number of Dofs.');

            obj.weights = cellfun(@(w) w(:), weights, 'Un', 0);
        end

        function obj = setWeight(obj, index, weight)
            % SETWEIGHT Set a specific weight.
            %
            %   obj = setWeight(obj, index, weight) sets the weight at the
            %   specified index.
            %
            % Inputs:
            %   obj - The SumSpace object
            %   index - Index of the factor (positive integer)
            %   weight - Weight vector (1 x nDofs vector)
            %
            % Outputs:
            %   obj - The SumSpace object
            
            tf = length(weight) == obj.getFactor(index).nDofs;
            core.except.assert(tf, 'InvalidInput', ...
                'Weight vector must match the number of Dofs.');

            obj.weights{index} = weight(:);
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the total number of degrees of freedom.

            n = sum(cellfun(@(F) F.nDofs, obj.factors));
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            if obj.nFactors > 0
                n = obj.factors{1}.nDims;
            else
                n = 0;
            end
        end

        function Y = evaluate(obj, X, varargin)
            % EVALUATE Evaluate functions in the direct sum space.
            %
            %   Y = evaluate(obj, X, C1, C2, ..., Ck) evaluates functions
            %   with coefficient matrices @a C1, @a C2, ..., @a Ck for each
            %   space at evaluation points @a X.
            %
            % Inputs:
            %   obj - The SumSpace object
            %   X - Evaluation points (nd x np matrix)
            %   varargin - Coefficients for each space:
            %<     Ci - Coefficients for space i (ni x nc matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points (nc x np matrix)

            X = reshape(X, size(X, 1), []);
            np = size(X, 2);
            m = obj.nFactors;

            C = varargin;
            nc = [];
            for i = 1:m
                C{i} = reshape(C{i}, size(C{i}, 1), []);
                if isempty(nc), nc = size(C{i}, 2); end
            end

            Y = zeros(nc, np);
            for i = 1:m
                D =  obj.weights{i} .* C{i};
                Z = obj.factors{i}.evaluate(X, D); % nc x np
                Y = Y + Z;
            end
        end

        function R = component(obj, i, Y, X, varargin)
            % COMPONENT Evaluate the i-th component.
            %
            %   R = component(obj, i, Y, X, C1, C2, ..., Ck) evaluates the
            %   @a i-th component with coefficient matrices @a C1, @a C2,
            %   ..., @a Ck for each space at evaluation points @a X.
            %
            % Inputs:
            %   obj - The SumSpace object
            %   Y - Target function values (nc x np matrix)
            %   X - Evaluation points (nd x np matrix)
            %   varargin - Coefficients for each space:
            %<     Ci - Coefficients for space i (ni x nc matrix)
            %
            % Outputs:
            %   Y - Function values at evaluation points (nc x np matrix)
            
            m = obj.nFactors;
            C = varargin;
            R = Y;
            for k = setdiff(1:m, i)
                D = obj.weights{k} .* C{k};
                Z = obj.factors{k}.evaluate(X, D); % nc x np
                R = R - Z;
            end
            R = R ./ obj.weights{i};
        end
    end
end