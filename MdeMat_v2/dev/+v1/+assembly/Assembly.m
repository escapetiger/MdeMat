classdef Assembly < handle
    % ASSEMBLY Base class for assembly procedure.
    %
    %   Assembly provides the foundation for assembling finite element
    %   operators into global operators. This abstract base class manages
    %   finite element spaces and operators, and provides utilities for
    %   sparse matrix construction.
    %
    %   The class handles the relationship between finite elements, meshes,
    %   and operators, providing common properties and methods needed for
    %   matrix assembly operations.
    %
    % See also:
    %   fem.assembly.ConstantAssembly, fem.assembly.GridAssembly

    properties
        fe   % Finite element object
        mesh % Mesh object  
        op   % Finite element operator object
    end

    properties (Dependent)
        nTotalElements % Total number of elements
        nLocalDofs     % Number of degrees of freedom per element
        nGlobalDofs    % Total number of degrees of freedom
        nDims          % Number of dimensions
    end

    methods
        function obj = Assembly(fe, mesh, op)
            % ASSEMBLY Constructor for Assembly.
            %
            %   obj = Assembly(fe, mesh, op) creates an assembler object
            %   with the specified finite element, mesh, and operator.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Mesh object
            %   op - Finite element operator object
            %
            % Outputs:
            %   obj - Constructed Assembly object

            obj.fe = fe;
            obj.mesh = mesh;
            obj.op = op;
        end

        function n = get.nTotalElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.

            n = obj.mesh.nTotalElements;
        end

        function n = get.nLocalDofs(obj)
            % GET.NLOCALDOFS Get the number of degrees of freedom per
            % element.

            n = obj.fe.nDofs;
        end

        function n = get.nGlobalDofs(obj)
            % GET.NGLOBALDOFS Get the total number of degrees of freedom.

            n = obj.nLocalDofs * obj.nTotalElements;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of dimensions.

            n = obj.fe.nDims;
        end
    end

    methods (Static)
        function A = sparseFromTriplets(T, m, n)
            % SPARSEFROMTRIPLETS Create sparse matrix from triplet format.
            %
            %   A = sparseFromTriplets(T, m, n) converts a triplet list 
            %   (row, column, value) into a sparse matrix, automatically
            %   filtering out invalid indices (zeros or negative values).
            %
            % Inputs:
            %   T - Triplet matrix with columns [rowIdx, colIdx, value]
            %   m - Number of rows in output sparse matrix (positive integer)
            %   n - Number of columns in output sparse matrix (positive integer)
            %
            % Outputs:
            %   A - Sparse matrix of size m×n

            k = find((T(:, 1) > 0) & (T(:, 2) > 0));
            A = sparse(T(k, 1), T(k, 2), T(k, 3), m, n);
        end
    end
end