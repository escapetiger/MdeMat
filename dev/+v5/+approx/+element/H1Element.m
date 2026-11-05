classdef H1Element < approx.element.Element
    % H1ELEMENT H1 Galerkin element.
    %
    %   H1Element implements a H1 element function. This class extends the
    %   base Element class by adding volume integration capabilities
    %   through precomputed basis values and first-order derivatives at
    %   integration points.
    %
    % Examples:
    %   % Create H1 element with quadrature
    %   geometry = core.geometry.Orthotope.unit(2);
    %   basis = approx.basis.LagrangeBasis(2);
    %   projector = approx.project.NodalProjector(basis);
    %   integrator = approx.integrate.GaussLegendre(3);
    %   element = H1Element(geometry, projector, integrator);
    %
    %   % Access pre-computed data
    %   values = element.volume.values;
    %   derivatives = element.volume.derivatives;
    %   weights = element.volume.weights;
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'h1' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
    end

    methods
        function obj = H1Element(geometry, projector, integrator)
            % H1ELEMENT Constructor for H1Element.
            %
            %   obj = H1Element(geometry, projector, integrator) creates a
            %   H1 element with the specified geometry, projector, and
            %   integrator for volume computations.
            %
            % Inputs:
            %   geometry - Geometry object
            %   projector - GalerkinProjector object
            %   integrator - Integrator object
            %
            % Outputs:
            %   obj - Constructed H1Element object

            obj@approx.element.Element(geometry, projector);
            obj.volume = approx.integrate.IntegrationData(integrator, 1, projector.basis);
        end
    end
end