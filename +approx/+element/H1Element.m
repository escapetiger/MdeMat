classdef H1Element < approx.element.L2Element
    % H1ELEMENT \f$H^1\f$ Sobolev space element.
    %
    %   H1Element implements an \f$H^1\f$ element function suitable for
    %   problems requiring both function values and first derivatives.
    %   This class extends the L2Element class by adding first-order
    %   derivative computation capabilities at integration points.
    %
    % See also:
    %   approx.integrate.Integrator, approx.element.L2Element

    methods
        function obj = H1Element(geometry, approximator, integrator)
            % H1ELEMENT Constructor for H1Element.
            %
            %   obj = H1Element(geometry, approximator, integrator) creates an
            %   H1 element with the specified @a geometry, @a approximator, and
            %   @a integrator for volume computations.

            arguments
                geometry core.geometry.Geometry
                approximator approx.linear.LinearApproximator
                integrator approx.integrate.Integrator
            end

            obj@approx.element.L2Element(geometry, approximator, integrator);
            obj.Volume = obj.Volume.setPartials(approximator.Basis, 1);
        end
    end
end