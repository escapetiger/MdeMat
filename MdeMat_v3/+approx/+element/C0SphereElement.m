classdef C0SphereElement < approx.element.C0Element
    % C0SPHEREELEMENT

    properties
        reduction % Reduction type
    end

    methods
        function obj = C0SphereElement(nDims, basisOrder, basisType, nIntPoints, reduction, radius)
            % SETC0SPHEREELEMENT Create a \f$C^0\f$ element on the sphere.
            %
            %   obj = C0SphereElement(nDims, basisOrder, basisType,
            %   nIntPoints, reduction, radius) creates a \f$C^0\f$ element
            %   on the sphere with specified basis configuration.
            %
            % Inputs:
            %   obj - The C0SphereElement object
            %   nDims - Number of dimensions (integer)
            %   basisOrder - Basis order (integer)
            %   basisType - Basis type (char)
            %   nIntPoints - Number of integration points (integer)
            %   reduction - Reduction type (1 or 2, default: 1)
            %   radius - Spere radius (positive scalar, default: 1)
            %
            % Outputs:
            %   obj - The C0SphereElement object
            %
            % Notes:
            %   Reduction parameter affects basis choice:
            %     - nDims=1, reduction=1: Legendre basis on {-1, 1}
            %     - nDims=1, reduction=2: Legendre basis on [-r, r]
            %     - nDims=2, reduction=1: Fourier basis on circle
            %     - nDims=3: Spherical harmonic basis on sphere

            if nargin < 5, reduction = 1; end
            if nargin < 6, radius = 1; end

            if nDims == 1 && reduction == 1
                A = 1 / 2;
                G = core.geometry.Hypersphere(0);
                E = core.function.ConstantFunction(1, A);
                I = approx.integrate.HypersphereIntegrator(1, ...
                    'gauss_trapezoidal');
                I.setPoints(nIntPoints);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
                switch basisType
                    case 'orthogonal'
                        F = approx.basis.CompiledOrthogonalBasis( ...
                            basisOrder, 'legendre', 'canonical');
                        F.autoLoad();
                    case 'cardinal'
                        F = approx.basis.CardinalKernelBasis(I.nodes);
                    case 'rbf'
                        F = approx.basis.RbfKernelBasis(I.nodes);
                end
            elseif nDims == 1 && reduction == 2
                A = 1 / (2 * radius);
                G = core.geometry.Orthotope([-1, 1]);
                E = core.function.ConstantFunction(1, A);
                I = approx.integrate.OrthotopeIntegrator( ...
                    1, 'gauss_legendre');
                I.setPoints(nIntPoints, -radius, radius);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
                switch basisType
                    case 'orthogonal'
                        F = approx.basis.CompiledOrthogonalBasis( ...
                            basisOrder, 'legendre', 'canonical'); 
                        F.autoLoad();
                    case 'cardinal'
                        F = approx.basis.CardinalKernelBasis(I.nodes);
                    case 'rbf'
                        F = approx.basis.RbfKernelBasis(I.nodes);
                end
            elseif nDims == 2 && reduction == 1
                A = 1 / (2 * pi * radius);
                G = core.geometry.Hypersphere(1);
                E = core.function.ConstantFunction(1, A);
                I = approx.integrate.HypersphereIntegrator( ...
                    2, 'gauss_trapezoidal');
                I.setPoints(nIntPoints, zeros(1, 2), radius, 3);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
                switch basisType
                    case 'orthogonal'
                        F = approx.basis.CompiledOrthogonalBasis( ...
                            basisOrder, 'fourier', 'canonical');
                        F.autoLoad();
                    case 'cardinal'
                        F = approx.basis.CardinalKernelBasis(I.nodes);
                    case 'rbf'
                        F = approx.basis.RbfKernelBasis(I.nodes);
                end
            else
                A = 1 / (4 * pi * radius^2);
                G = core.geometry.Hypersphere(2);
                E = core.function.ConstantFunction(2, A);
                I = approx.integrate.HypersphereIntegrator( ...
                    3, 'gauss_trapezoidal');
                I.setPoints(n, zeros(1, 3), radius, 3);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
                switch basisType
                    case 'orthogonal'
                        F = approx.basis.CompiledOrthogonalBasis( ...
                            basisOrder, 'spherical_harmonic', 'canonical');
                        F.autoLoad();
                    case 'cardinal'
                        F = approx.basis.CardinalKernelBasis(I.nodes);
                    case 'rbf'
                        F = approx.basis.RbfKernelBasis(I.nodes);
                end
            end
            B = core.function.ProductFunction(F, E);
            P = approx.linear.ModalApproximator(B);
            P.setMass(I.nodes, I.weights);
            P.setDesign(I.nodes, I.weights);
            obj@approx.element.C0Element(G, P, I);
            obj.reduction = reduction;
        end
    end
end
