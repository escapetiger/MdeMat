classdef HermiteState < physics.SpatialState
    % HERMITESTATE

    properties
        nModes
        history
    end

    methods
        function obj = HermiteState(xDisc, nh)
            obj@physics.SpatialState(xDisc);
            obj.dofs = struct( ...
                'D', [], ... %< Hermite coefficients
                'omega', [] ... %< Potential
                );
            obj.coefs = struct( ...
                'EInf', [], ... %< Electrical field steady state
                'rhoInf', [] ... %< Density steady state
                );
            obj.history = struct( ...
                't', [], ...
                'mass', [], ...
                'momentum', [], ...
                'energy', [], ...
                'potential', []);
            obj.nModes = nh + 1;
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined spatial state.
            %
            %   newObj = refine(obj, nLevels) creates a refined spatial
            %   state.
            %
            % Inputs:
            %   obj - The SpatialState object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined SpatialState object

            newXDisc = obj.xDisc.refine(nLevels);
            newObj = physics.vlasov.HermiteState(newXDisc, obj.nModes-1);
        end

        function [F, V] = distribution(obj, xRef, T0, sqrtRhoInfFn)
            NV = 512;
            VL = 3;
            V = linspace(-VL, VL, NV+1);
            D = obj.dofs.D;
            D = obj.xDisc.evaluate([], xRef, D);
            H = obj.hermite(V, T0);
            sqrtRhoInf = obj.xDisc.functionEvaluate([], xRef, sqrtRhoInfFn);
            F = (sqrtRhoInf .* D) * H;
        end

        function H = hermite(obj, V, T0)
            X = V ./ sqrt(T0);
            H = zeros(obj.nModes, size(X, 2));
            H(1, :) = 1;
            H(2, :) = X;
            M = exp(-V.^2./2) / sqrt(2*pi);
            for k = 2:obj.nModes - 1
                H(k+1, :) = (X .* H(k, :) - H(k-1, :) * sqrt(k-1)) / sqrt(k);
            end
            H = H .* M;
        end
    
        function obj = addHistory(obj, t, T0, rhoInf)
            obj.history.t(end+1) = t;

            D0 = obj.dofs.D(:, 1);
            D1 = obj.dofs.D(:, 2);
            D2 = obj.dofs.D(:, 3);
            omega = obj.xDisc.evaluate([], [], obj.dofs.omega);
            omega = omega(:).';
            dOmega = obj.xDisc.differentiate([], [], obj.dofs.omega, 1);
            dOmega = dOmega(:).';
            psi = @(x) (T0 * (dOmega - omega ./ rhoInf(x)) ./ sqrt(rhoInf(x))).^2;

            n = obj.xDisc.element.volume.nPoints;
            w = obj.xDisc.element.volume.weights;
            xRef = obj.xDisc.element.volume.nodes;
            h = obj.xDisc.mesh.computeElementJacobianDeterminants();

            mass = obj.xDisc.evaluate([], xRef, D0);
            mass = w * reshape(mass, n, []);
            mass = sum(mass .* h(:).');
            obj.history.mass(end+1) = mass;

            momentum = obj.xDisc.evaluate([], xRef, D1);
            momentum = w * reshape(momentum, n, []);
            momentum = sum(momentum.* h(:).');
            obj.history.momentum(end+1) = momentum;

            energy = obj.xDisc.evaluate([], xRef, (sqrt(2)*D2+D0)/2);
            energy = w * reshape(energy, n, []);
            energy = sum(energy.* h(:).');
            obj.history.energy(end+1) = energy;

            potential = obj.xDisc.functionEvaluate([], xRef, psi);
            potential = w * reshape(potential, n, []);
            potential = sum(potential.* h(:).');
            obj.history.potential(end+1) = potential;
        end
    end
end
