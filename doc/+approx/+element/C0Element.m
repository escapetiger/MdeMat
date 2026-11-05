classdef C0Element < approx.element.Element
    % C0ELEMENT Element-level interior-continuous function.
    %
    %   C0Element implements an element-level function that is continuous
    %   across element interior and supports volume integration. This class
    %   extends the base Element class by adding volume integration
    %   capabilities through pre-computed function data at quadrature
    %   points.
    %
    %   The class manages function evaluations and derivatives at volume
    %   integration points, enabling efficient assembly of finite element
    %   matrices and vectors. It is designed for conforming finite element
    %   methods where continuity across element boundaries is enforced.
    %
    % Examples:
    %   % Create C0 element with quadrature
    %   geometry = core.geometry.Orthotope.unit(2);
    %   basis = approx.basis.LagrangeBasis(2);
    %   projector = approx.project.NodalProjector(basis);
    %   integrator = approx.integrate.GaussLegendre(3);
    %   element = C0Element(geometry, projector, integrator);
    %
    %   % Set up function data for mass matrix assembly
    %   element.setVolumeData(0);  % Function values only
    %
    %   % Set up function data for stiffness matrix assembly
    %   element.setVolumeData(1);  % Include first derivatives
    %
    %   % Access pre-computed data
    %   values = element.volumeData.values;
    %   weights = element.volumeData.weights;
    %
    % See also:
    %   approx.element.Element, approx.element.L2Element,
    %   approx.element.ElementFunction

    properties
        volumeData % Function data for volume integrals (ElementFunction)
    end

    methods
        function obj = C0Element(geometry, projector, integrator)
            % C0ELEMENT Constructor for C0Element.
            %
            %   obj = C0Element(geometry, projector, integrator) creates a
            %   continuous element with the specified geometry, projector,
            %   and integrator for volume computations.
            %
            % Inputs:
            %   geometry - Geometry object defining the element domain
            %   projector - GalerkinProjector object for basis functions
            %   integrator - Integrator object for volume quadrature
            %
            % Outputs:
            %   obj - Constructed C0Element object

            obj@approx.element.Element(geometry, projector);
            obj.volumeData = approx.element.ElementFunction(integrator);
        end

        function obj = setVolumeData(obj, k)
            % SETVOLUMEDATA Set function data for volume integrals.
            %
            %   obj = setVolumeData(obj, k) computes and stores function
            %   values and derivatives up to order @a k at integration
            %   points for volume integral computations. This method
            %   evaluates the basis functions and their derivatives at
            %   quadrature points, enabling efficient matrix assembly.
            %
            % Inputs:
            %   obj - The C0Element object
            %   k - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   obj - The C0Element object

            D = obj.volumeData;
            X = D.integrator.nodes(1:obj.nDims, :);
            D.setValues(obj.projector.basis, X);
            if k > 0, D.setDerivatives(obj.projector.basis, X, k); end
        end
    end
end