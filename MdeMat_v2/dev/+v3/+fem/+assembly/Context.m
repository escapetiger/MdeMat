classdef Context < handle
    % ASSEMBLY Base class for assembly contexts.
    %
    %   Context encapsulates finite element, mesh and finite element
    %   operator to provide common properties and methods in assembly
    %   procedures.
    %
    % See Also:
    %   fem.assembly.Context

    properties
        fe % Finite element object
        mesh % Mesh object
        feOp % Finite element operator object
    end

    properties (Dependent)
        nTotalElements % Total number of elements
        nLocalDofs % Number of degrees of freedom per element
        nGlobalDofs % Total number of degrees of freedom
        nDims % Number of dimensions
    end

    methods
        function obj = Context(fe, mesh, feOp)
            % CONTEXT Constructor for Context.
            %
            %   obj = Context(fe, mesh, feOp) creates an assembly context
            %   with the specified finite element, mesh, and operator.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Mesh object
            %   feOp - Finite element operator object
            %
            % Outputs:
            %   obj - Constructed Context object

            obj.fe = fe;
            obj.mesh = mesh;
            obj.feOp = feOp;
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
end