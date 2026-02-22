classdef FiniteElementOperator < handle
    % FINITEELEMENTOPERATOR Base class for finite element operators.
    %
    %   FiniteElementOperator provides a foundation for operators that act
    %   on finite element spaces, such as gradient, divergence, curl, and
    %   other differential operators. This abstract base class maintains
    %   a reference to the underlying finite element and provides common
    %   functionality for all derived operator classes.
    %
    %   Finite element operators encapsulate the discrete forms of
    %   differential operators, enabling efficient computation of terms
    %   in weak formulations of partial differential equations. They
    %   bridge the gap between continuous differential operators and
    %   their discrete matrix representations.
    %
    % See also:
    %   fem.element.C0FiniteElementOperator, fem.element.L2FiniteElementOperator,
    %   fem.element.SemiLagrangianFiniteElementOperator

    properties
        fe % Finite element object on which the operator acts
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom (dependent property)
    end

    methods
        function obj = FiniteElementOperator(fe)
            % FINITEELEMENTOPERATOR Constructor for FiniteElementOperator.
            %
            %   obj = FiniteElementOperator(fe) creates a finite element
            %   operator associated with the specified finite element.
            %   This constructor is called by concrete subclasses to
            %   initialize the base class functionality.
            %
            % Inputs:
            %   fe - Finite element object (FiniteElement)
            %
            % Outputs:
            %   obj - Constructed FiniteElementOperator object

            obj.fe = fe;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.

            n = obj.fe.nDofs;
        end
    end
end