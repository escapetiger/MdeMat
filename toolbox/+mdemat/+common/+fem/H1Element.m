classdef H1Element < mdemat.common.fem.L2Element
    % H1ELEMENT \f$H^1\f$ Sobolev space element.
    %
    %   H1Element implements an \f$H^1\f$ element function suitable for
    %   problems requiring both function values and first derivatives.
    %   This class extends the L2Element class by adding first-order
    %   derivative computation capabilities at integration points.
    %
    % See also:
    %   mdemat.common.quadrature.Integrator, mdemat.common.fem.L2Element

    methods
        function obj = H1Element(geometry, approximator, integrator)
            % H1ELEMENT Constructor for H1Element.
            %
            %   obj = H1Element(geometry, approximator, integrator) creates an
            %   H1 element with the specified @a geometry, @a approximator, and
            %   @a integrator for volume computations.

            arguments
                geometry mdemat.data.geometry.Geometry
                approximator mdemat.common.approximation.LinearApproximator
                integrator mdemat.common.quadrature.Integrator
            end

            obj@mdemat.common.fem.L2Element(geometry, approximator, integrator);
            obj.Volume = obj.Volume.setPartials(approximator.Basis, 1);
        end
    end
end