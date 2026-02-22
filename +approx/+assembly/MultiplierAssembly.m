classdef MultiplierAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTMULTIPLIERASSEMBLY Assembly for finite element
    % multipliers.
    %
    %   MultiplierAssembly provides functionality for
    %   assembling matrices from multipliers in finite element
    %   formulations. This class handles both constant coefficient and
    %   spatially varying coefficient multiplication operators.
    %
    % See also:
    %   approx.assembly.FiniteElementAssembly

    methods
        function G = assembleMatrix(obj, field, options)
            % ASSEMBLEMATRIX Assemble matrix with a fixed field.
            %
            %   G = assembleMatrix(obj, field) assembles the global matrix
            %   for the multiplication operator with field coefficients @a field.
            %
            %   G = assembleMatrix(obj, field, options) passes additional
            %   arguments through @a options.args to function evaluation.

            arguments
                obj approx.assembly.MultiplierAssembly
                field{mustBeA(field, {'numeric', 'function_handle'})}
                options.args = {}
            end

            if isa(field, 'function_handle')
                G = obj.assemblePointwiseMatrix(field, options.args{:});
            else
                G = obj.assembleConstantMatrix(field);
            end
        end
    end

    methods (Access = protected)
        function G = assembleConstantMatrix(obj, field)
            % ASSEMBLECONSTANTMATRIX Assemble matrix with constant field.
            %
            %   G = assembleConstantMatrix(obj, field) assembles the global
            %   matrix for the multiplication operator with constant field
            %   coefficients @a field using a diagonal matrix representation.

            arguments
                obj approx.assembly.MultiplierAssembly
                field {mustBeNumeric}
            end

            ng = obj.Space.NGlobalDofs;
            nl = obj.Space.NLocalDofs;

            if isscalar(field)
                field = repmat(field, obj.Space.NMeshElements, 1);
            end

            G = spdiags(kron(field(:), ones(nl, 1)), 0, ng, ng);
        end

        function G = assemblePointwiseMatrix(obj, field, varargin)
            % ASSEMBLEPOINTWISEMATRIX Assemble matrix with spatially
            % varying field.
            %
            %   G = assemblePointwiseMatrix(obj, field, varargin) assembles
            %   the global matrix for the multiplication operator with
            %   spatially varying field coefficients @a field using volume
            %   integration with weighted mass matrices.

            arguments
                obj approx.assembly.MultiplierAssembly
                field {mustBeA(field, 'function_handle')}
            end
            arguments (Repeating)
                varargin
            end

            nq = obj.Space.Element.Volume.NPoints;
            ne = obj.Space.NMeshElements;
            xRef = obj.Space.Element.Volume.Nodes;
            A = obj.Space.feval(field, xRef, args = varargin);
            A = reshape(A, nq, []);
            W0 = obj.Space.Element.Volume.Weights; % (1, nq)
            W0 = W0(:) .* A; % (nq, ne)
            W = zeros(nq, nq, ne);
            [i, j] = ndgrid(1:nq, 1:ne);
            I = sub2ind([nq, nq, ne], i, i, j);
            W(I) = W0(:);
            V = obj.Space.Element.Volume.Values; % (nl, nq)
            U = obj.Space.Element.Volume.Values; % (nl, nq)
            LA = pagemtimes(pagemtimes(V, W), U.');
            LA = obj.Space.Element.Approximator.project(LA);
            LA = reshape(LA, size(V, 1), size(U, 1), []);
            G = obj.createMatrix(LA);
        end
    end
end
