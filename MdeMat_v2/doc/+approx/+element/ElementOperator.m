classdef ElementOperator < handle
    % ELEMENTOPERATOR Base class for element operators.
    %
    %   ElementOperator encapsulates the discrete forms of differential
    %   operators, bridging the gap between continuous differential
    %   operators and their discrete representations on finite elements.
    %
    %   The class serves as a common interface for operators that act on
    %   element-level data, including gradients, divergences, and other
    %   differential operators. Concrete subclasses implement specific
    %   operator types for different element conformity requirements
    %   (C0, L2) and discretization approaches (standard, semi-Lagrangian).
    %
    % See also:
    %   approx.element.C0ElementOperator, 
    %   approx.element.L2ElementOperator,
    %   approx.element.SemiLagrangianElementOperator

    properties
        element % Element object on which the operator acts
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom (integer)
    end

    methods
        function obj = ElementOperator(element)
            % ELEMENTOPERATOR Constructor for ElementOperator.
            %
            %   obj = ElementOperator(element) creates an operator
            %   associated with the specified element. This constructor
            %   establishes the connection between the operator and its
            %   underlying element representation.
            %
            % Inputs:
            %   element - Element object on which the operator will act
            %
            % Outputs:
            %   obj - Constructed ElementOperator object

            obj.element = element;
        end

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.

            n = obj.element.nDofs;
        end
    end
end