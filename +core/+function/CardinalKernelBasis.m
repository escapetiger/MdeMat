classdef CardinalKernelBasis < core.function.Function
    % CARDINALKERNELBASIS Interpolation basis using cardinal kernel
    % methods.
    %
    %   CardinalKernelBasis represents a set of cardinal interpolation
    %   basis functions using barycentric coordinates for efficient
    %   evaluation. The basis functions are constructed from given node
    %   points and provide interpolation properties at those nodes.
    %
    % See Also:
    %   core.function.Function

    properties
        Nodes % Interpolation nodes (nd x nc)
        Weights % Barycentric weights (nc x 1)
    end

    methods
        function obj = CardinalKernelBasis(nodes)
            % CARDINALKERNELBASIS Constructor for CardinalKernelBasis.
            %
            %   obj = CardinalKernelBasis(nodes) creates a
            %   CardinalKernelBasis object from the specified node points
            %   and computes barycentric weights for efficient evaluation.
            %   Duplicate nodes are automatically removed to ensure
            %   well-posed interpolation.

            arguments
                nodes {mustBeNumeric}
            end

            nDims = size(nodes, 1);
            nCodims = size(nodes, 2);

            obj@core.function.Function(nDims=nDims, nCodims=nCodims);
            obj.setNodes(nodes);
        end

        function obj = setNodes(obj, nodes)
            % SETNODES Set the interpolation nodes.
            %
            %   obj = setNodes(nodes) sets interpolation nodes and update
            %   barycentric weights.

            arguments
                obj
                nodes {mustBeNumeric}
            end

            obj.Nodes = nodes;
            obj.updateWeights(nodes);
        end
    end

    methods (Access = protected)
        function obj = updateWeights(obj, nodes)
            % SETWEIGHTS Compute barycentric weights for interpolation.
            %
            %   obj = updateWeights(obj, nodes) computes the barycentric
            %   weights used in the barycentric cardinal interpolation
            %   formula for efficient evaluation. Weights are computed
            %   as reciprocals of products of inter-node distances.

            X1 = reshape(nodes, [obj.NDims, obj.NCodims, 1]);
            X2 = reshape(nodes, [obj.NDims, 1, obj.NCodims]);
            DX = vecnorm(X1-X2, 2, 1); %< [1, nc, nc]
            DX = reshape(DX, [obj.NCodims, obj.NCodims]); %< [nc, nc]
            DX(DX == 0) = 1;
            obj.Weights = 1 ./ prod(DX, 2); %< [nc, 1]
        end

        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of cardinal basis evaluation.
            %
            %   Y = evalImpl(obj, X) evals the cardinal basis functions at
            %   the specified points using the barycentric formula for
            %   efficient computation. The method handles exact node
            %   matches to avoid division by zero.

            nd = obj.NDims;
            nc = obj.NCodims;
            np = size(X, 2);
            X0 = reshape(obj.Nodes, [nd, nc, 1]); %< [nd, nc, 1]
            X0 = repmat(X0, [1, 1, np]); %< [nd, nc, np]
            X = reshape(X, [nd, 1, np]); %< [nd, 1, np]
            X = repmat(X, [1, nc, 1]); %< [nd, nc, np]
            DX = vecnorm(X-X0, 2, 1); %< [nc, np]
            DX = reshape(DX, [nc, np]); %< [nc, np]

            %< apply barycentric formula
            W = repmat(obj.Weights, [1, np]); %< [nc, np]
            P = W ./ DX; %< [nc, np]
            Y = P ./ sum(P, 1); %< normalized, [nc, np]
            Y(DX == 0) = 1;
        end
    end
end