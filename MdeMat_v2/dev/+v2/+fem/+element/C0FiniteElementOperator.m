classdef C0FiniteElementOperator < fem.element.FiniteElementOperator
    % C0FINITEELEMENTOPERATOR Operator for continuous finite elements.
    %
    %   C0FiniteElementOperator provides differential operators such as
    %   gradient for C0 finite elements. It manages bilinear form data
    %   for volume integrals and constructs discrete operator matrices
    %   used in conforming finite element methods.
    %
    %   The class handles the assembly of discrete differential operators
    %   by pre-computing bilinear forms involving function derivatives
    %   at quadrature points. This enables efficient construction of
    %   system matrices for elliptic and parabolic problems.
    %
    %   For C0 elements, operators only involve volume integrals since
    %   function continuity across element boundaries eliminates the
    %   need for explicit flux terms in the weak formulation.
    %
    % Examples:
    %   % Create C0 element and its operator
    %   element = fem.element.C0FiniteElement(geometry, projector, integrator);
    %   element.setVolumeData(1);  % Include first derivatives
    %   operator = C0FiniteElementOperator(element);
    %   
    %   % Set up operator data and access gradient
    %   operator.setVolumeData();
    %   gradOp = operator.gradient;
    %   
    %   % Use gradient components for stiffness matrix assembly
    %   gradX = gradOp.volumeData(1).matrix;  % ∂/∂x component
    %   gradY = gradOp.volumeData(2).matrix;  % ∂/∂y component
    %
    % See also:
    %   fem.element.FiniteElementOperator, fem.element.L2FiniteElementOperator,
    %   fem.element.C0FiniteElement
    
    properties
        volumeData % Array of bilinear form data for volume integrals (BilinearFormData array)
    end

    properties (Dependent)
        gradient % Gradient operator matrices (dependent property)
    end

    methods
        function obj = C0FiniteElementOperator(fe)
            % C0FINITEELEMENTOPERATOR Constructor for
            % C0FiniteElementOperator.
            %
            %   obj = C0FiniteElementOperator(fe) creates a continuous
            %   finite element operator for the specified C0 finite
            %   element. The operator provides access to discrete
            %   differential operators through pre-computed bilinear forms.
            %
            % Inputs:
            %   fe - C0 finite element object (C0FiniteElement)
            %
            % Outputs:
            %   obj - Constructed C0FiniteElementOperator object
            
            cls = 'fem.element.C0FiniteElement';
            core.except.assert(isa(fe, cls), ...
                'InvalidInput', 'Input must be a C0 FiniteElement.');

            obj@fem.element.FiniteElementOperator(fe);
        end

        function G = get.gradient(obj)
            % GET.GRADIENT Get the gradient operator.
            
            G = obj.buildGradient();
        end

        function obj = setVolumeData(obj)
            % SETVOLUMEDATA Set up bilinear form data for volume integrals.
            %
            %   obj = setVolumeData(obj) creates bilinear form matrices for
            %   each derivative order available in the finite element's
            %   volume data. These matrices represent discrete differential
            %   operators computed using numerical integration.
            %
            % Inputs:
            %   obj - The C0FiniteElementOperator object
            %
            % Outputs:
            %   obj - The C0FiniteElementOperator object
            
            D = obj.fe.volumeData;
            m = D.nDerivatives;
            
            %< Create bilinear forms for each derivative order
            obj.volumeData = arrayfun( ...
                @(i) fem.element.BilinearFormData(D, D, i, 0), 1:m);
        end
    end

    methods (Access = protected)
        function G = buildGradient(obj)
            % BUILDGRADIENT Build gradient operator matrices.
            %
            %   G = buildGradient(obj) extracts the gradient components
            %   from the volume data corresponding to first-order derivatives.
            %   The gradient operator maps scalar functions to vector fields.
            %
            % Inputs:
            %   obj - The C0FiniteElementOperator object
            %
            % Outputs:
            %   G - Structure with volume gradient data
            %       .volumeData - Volume gradient operators for each dimension
            
            nd = obj.fe.nDims;
            G = struct('volumeData', obj.volumeData(1:nd));
        end
    end
end