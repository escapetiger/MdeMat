classdef L2Element < approx.element.Element
    % L2ELEMENT Continuous element for \f$L^2\f$ function spaces.
    %
    %   L2Element implements a continuous element function for
    %   approximating functions in \f$L^2\f$ spaces. This class extends
    %   the base Element class by adding volume integration capabilities
    %   through precomputed basis function values at integration points.
    %
    %   The \f$L^2\f$ element serves as the base class for elements that
    %   require volume integration, making it suitable for various finite
    %   element methods and spectral methods.
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties
        Volume % Volume integrator (Integrator)
    end

    methods
        function obj = L2Element(geometry, approximator, integrator)
            % L2ELEMENT Constructor for L2Element.
            %
            %   obj = L2Element(geometry, approximator, integrator)
            %   creates an L2 element with the specified @a geometry,
            %   @a approximator, and @a integrator for volume computations.

            arguments
                geometry core.geometry.Geometry
                approximator approx.linear.LinearApproximator
                integrator approx.integrate.Integrator
            end

            obj@approx.element.Element(geometry, approximator);
            obj.Volume = integrator.setValues(approximator.Basis);
        end
    end
end