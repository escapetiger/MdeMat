classdef C0Element < approx.element.Element
    % C0ELEMENT Continuous element for \f$C^0\f$ function spaces.
    %
    %   C0Element implements a continuous element function for
    %   approximating functions in \f$C^0\f$ spaces (continuous functions).
    %   This class extends the base Element class by adding volume
    %   integration capabilities through precomputed basis function values
    %   at integration points.
    %
    %   The \f$C^0\f$ element ensures continuity of the approximated
    %   function across element boundaries, making it suitable for problems
    %   requiring global continuity such as displacement-based finite
    %   element methods or spectral methods.
    %
    % See also:
    %   approx.integrate.IntegrationData, approx.element.Element

    properties (Constant)
        TYPE = 'c0' % Element type identifier
    end

    properties
        volume % Volume integration data (IntegrationData)
    end

    methods
        function obj = C0Element(geometry, approximator, integrator)
            % C0ELEMENT Constructor for C0Element.
            %
            %   obj = C0Element(geometry, approximator, integrator)
            %   creates a continuous element with the specified geometry,
            %   approximator, and integrator for volume computations.
            %
            % Inputs:
            %   geometry - Geometry object
            %   approximator - LinearApproximator object
            %   integrator - Integrator object
            %
            % Outputs:
            %   obj - Constructed C0Element object

            obj@approx.element.Element(geometry, approximator);
            obj.volume = approx.integrate.IntegrationData(integrator, 0, approximator.basis);
        end
    end
end