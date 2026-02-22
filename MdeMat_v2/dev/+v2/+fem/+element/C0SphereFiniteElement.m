classdef C0SphereFiniteElement < fem.element.C0FiniteElement
    % C0SPHEREFINITEELEMENT

    methods
        function obj = C0SphereFiniteElement(varargin)
            % C0SPHEREFINITEELEMENT
            p = inputParser;
            addParameter(p, 'nDims', [], @(x) x >= 1 && x <= 3);
            addParameter(p, 'reduction', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'order', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'radius', 1, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'nPoints', [], @(x) isnumeric(x));
            parse(p, varargin{:});

            r = p.Results.radius;
            if p.Results.nDims == 1 && p.Results.reduction == 1
                K = core.geometry.Hypersphere(0);
                F = approx.basis.OrthogonalBasisFunction( ...
                    p.Results.order, 'legendre', 'canonical');
                E = core.function.ConstantFunction(1, 1/2);
                I = approx.integrate.HypersphereIntegrator(1, 'gauss_trapezoidal');
                I.setPoints(p.Results.nPoints, false);
            elseif p.Results.nDims == 1 && p.Results.reduction == 2
                K = core.geometry.Orthotope([-1, 1]);
                F = approx.basis.OrthogonalBasisFunction( ...
                    p.Results.order, 'legendre', 'canonical');
                E = core.function.ConstantFunction(1, 1/(2 * r));
                I = approx.integrate.OrthotopeIntegrator(1, 'gauss_legendre');
                I.setPoints(p.Results.nPoints, false, -1, 1);
            elseif p.Results.nDims == 2 && p.Results.reduction == 1
                K = core.geometry.Hypersphere(1);
                F = approx.basis.OrthogonalBasisFunction( ...
                    p.Results.order, 'fourier', 'canonical');
                E = core.function.ConstantFunction(1, 1/(2 * pi * r));
                I = approx.integrate.HypersphereIntegrator(2, 'gauss_trapezoidal');
                I.setPoints(p.Results.nPoints, false, zeros(1, 2), r, 3);
            else
                K = core.geometry.Hypersphere(2);
                F = approx.basis.OrthogonalBasisFunction( ...
                    p.Results.order, 'spherical_harmonic', 'canonical');
                E = core.function.ConstantFunction(2, 1/(4 * pi * r^2));
                I = approx.integrate.HypersphereIntegrator(3, 'gauss_trapezoidal');
                I.setPoints(p.Results.nPoints, false, zeros(1, 3), r, 3);
            end
            F.autoLoad();
            I.addWeightFunction(@(v) 1./E.evaluate(v));
            B = core.function.ProductFunction(F, E);
            P = approx.project.ModalProjector(B);
            P.setMass(I.nodes, I.weights);
            obj@fem.element.C0FiniteElement(K, P, I);
        end
    end
end
