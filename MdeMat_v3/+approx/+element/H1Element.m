classdef H1Element < approx.element.Element
    % H1ELEMENT \f$H^1\f$ Sobolev space element.
    %
    %   H1Element implements an \f$H^1\f$ element function suitable for
    %   problems requiring both function values and first derivatives.
    %   This class extends the base Element class by adding volume 
    %   integration capabilities through precomputed basis function values
    %   and first-order derivatives at integration points.
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
        function obj = H1Element(geometry, approximator, integrator)
            % H1ELEMENT Constructor for H1Element.
            %
            %   obj = H1Element(geometry, approximator, integrator) creates a
            %   H1 element with the specified geometry, approximator, and
            %   integrator for volume computations.
            %
            % Inputs:
            %   geometry - Geometry object
            %   approximator - approximator object
            %   integrator - Integrator object
            %
            % Outputs:
            %   obj - Constructed H1Element object

            obj@approx.element.Element(geometry, approximator);
            obj.volume = approx.integrate.IntegrationData(integrator, 1, approximator.basis);
        end
    end
end