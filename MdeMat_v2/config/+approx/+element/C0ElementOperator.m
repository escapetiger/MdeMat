classdef C0ElementOperator < approx.element.ElementOperator
    % C0ELEMENTOPERATOR Operator for continuous finite elements.
    %
    %   C0ElementOperator provides differential operators such as gradient
    %   for C0 elements. It manages bilinear form data for volume integrals
    %   and constructs discrete operator matrices used in conforming finite
    %   element methods.
    %
    %   The class computes discrete representations of differential
    %   operators by evaluating weighted inner products between basis
    %   function derivatives. These operators are fundamental for
    %   assembling system matrices in finite element discretizations of
    %   partial differential equations.
    %
    % Examples:
    %   % Create C0 element and its operator
    %   element = approx.element.C0Element(geometry, projector, integrator);
    %   element.setVolumeData(1);  % Include first derivatives
    %   operator = C0ElementOperator(element);
    %   
    %   % Set up operator data and access gradient
    %   operator.setVolumeData();
    %   gradOp = operator.gradient;
    %   
    %   % Use gradient components for stiffness matrix assembly
    %   gradX = gradOp.volumeData(1).matrix;  % x-derivative component
    %   gradY = gradOp.volumeData(2).matrix;  % y-derivative component
    %
    %   % Access complete gradient structure
    %   nDims = length(gradOp.volumeData);
    %   stiffnessContrib = sum(arrayfun(@(i) gradOp.volumeData(i).matrix, 1:nDims), 3);
    %
    % See also:
    %   approx.element.ElementOperator, 
    %   approx.element.L2ElementOperator,
    %   approx.element.C0Element
    
    properties
        volumeData % Array of bilinear form data for volume integrals (ElementBilinearForm array)
    end

    properties (Dependent)
        gradient % Gradient operator matrices (struct with volumeData field)
    end

    methods
        function obj = C0ElementOperator(element)
            % C0ELEMENTOPERATOR Constructor for C0ElementOperator.
            %
            %   obj = C0ElementOperator(element) creates an element
            %   operator for the specified C0 element. The element must
            %   have derivative data available for operator construction.
            %
            % Inputs:
            %   element - C0 element object (C0Element)
            %
            % Outputs:
            %   obj - Constructed C0ElementOperator object
            
%             cls = 'approx.element.C0Element';
%             core.except.assert(isa(element, cls), ...
%                 'InvalidInput', 'Input must be a C0 Element.');

            obj@approx.element.ElementOperator(element);
        end

        function G = get.gradient(obj)
            % GET.GRADIENT Get the gradient operator.
            
            G = obj.extractGradient();
        end

        function obj = setVolumeData(obj)
            % SETVOLUMEDATA Set up bilinear form data for volume integrals.
            %
            %   obj = setVolumeData(obj) creates discrete bilinear forms
            %   representing volume integrals for each derivative order.
            %   These forms are used to construct differential operators
            %   such as gradient and divergence for finite element
            %   computations.
            %
            % Inputs:
            %   obj - The C0ElementOperator object
            %
            % Outputs:
            %   obj - The C0ElementOperator object
            
            D = obj.element.volumeData;
            m = D.nDerivatives;
            f = @(i) approx.element.ElementBilinearForm(D, D, i, 0);
            obj.volumeData = arrayfun(f, 1:m);
        end
    end

    methods (Access = protected)
        function G = extractGradient(obj)
            % EXTRACTGRADIENT Extract gradient matrices.
            %
            %   G = extractGradient(obj) extracts the gradient data from
            %   the volume data corresponding to first-order partial
            %   derivatives. The number of components equals the spatial
            %   dimension of the element.
            %
            % Inputs:
            %   obj - The C0ElementOperator object
            %
            % Outputs:
            %   G - Gradient data structure with volumeData field
            
            G = struct('volumeData', obj.volumeData(1:obj.element.nDims));
        end
    end
end