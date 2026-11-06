classdef TestHermiteForce < matlab.unittest.TestCase
    %TESTHERMITEFORCE Compare mode recursions against global block solves.
    methods (Test)
        function standardStagesSolveTheirFullTriangularSystem(testCase)
            [xSpace, vSpace, coefficients, multiplier] = ...
                TestHermiteForce.fixture();
            state = mdemat.application.plasma.HermiteState(xSpace, vSpace);
            state.Dofs.E = ones(xSpace.NGlobalDofs, 1);
            scheme = HermiteForceProbe(config=struct('epsilon',0.37,'T0',1.7));
            scheme.Multiplier = FixedMultiplier(multiplier);
            methods = {'be', 'midpoint', 'sdirk2'};
            dt = 0.23;
            for i = 1:numel(methods)
                state.Dofs.D = coefficients;
                scheme.forceStep(state, dt, methods{i});
                expected = TestHermiteForce.blockSolve(coefficients, ...
                    multiplier/(0.37*sqrt(1.7)), zeros(size(coefficients)), ...
                    dt, methods{i});
                testCase.verifyEqual(state.Dofs.D, expected, ...
                    'AbsTol', 5e-13);
            end
        end

        function equilibriumStagesSolveTheirAffineSystem(testCase)
            [xSpace, vSpace, coefficients, multiplier] = ...
                TestHermiteForce.fixture();
            state = mdemat.application.plasma.ApHermiteState(xSpace, vSpace);
            ng = xSpace.NGlobalDofs;
            state.Dofs.Q = linspace(0.4, 1.2, ng).';
            state.Coefs.sqrtRhoInf = linspace(0.8, 1.1, ng).';
            equilibrium = zeros(size(coefficients));
            equilibrium(:, 1) = state.Coefs.sqrtRhoInf;
            scheme = ApHermiteForceProbe(config=struct('epsilon',0.29));
            scheme.MultAmb = FixedMultiplier(multiplier);
            methods = {'be', 'sdirk2'};
            dt = 0.17;
            for i = 1:numel(methods)
                state.Dofs.D = coefficients;
                scheme.forceStep(state, dt, methods{i});
                expected = TestHermiteForce.blockSolve(coefficients, ...
                    multiplier/0.29, equilibrium, dt, methods{i});
                testCase.verifyEqual(state.Dofs.D, expected, ...
                    'AbsTol', 5e-13);
            end
            state.Dofs.D = equilibrium;
            scheme.forceStep(state, dt, 'be');
            testCase.verifyEqual(state.Dofs.D, equilibrium, 'AbsTol', 1e-14);
        end

        function equilibriumGaugeUsesPhysicalWeightedIntegral(testCase)
            mesh = mdemat.data.mesh.Grid({[0, 0.12, 0.6, 1]});
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 2);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            velocity = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1, 5));
            state = mdemat.application.plasma.ApHermiteState(space, velocity);
            scheme = ApHermiteForceProbe(config=struct( ...
                'sqrtRhoInf',@(x) 1+0.2*cos(2*pi*x(1,:))));
            scheme.MultAmb = mdemat.common.fem.MultiplierAssembly(space);
            scheme.SrcAmb = mdemat.common.fem.SourceAssembly(space);
            [constant, weights] = scheme.gaugeWeights(state);
            coefficients = sin((1:space.NGlobalDofs).');
            nodes = element.Volume.Nodes;
            physical = mesh.collocate(nodes);
            weighted = space.eval(nodes, coefficients) ...
                ./ scheme.Config.sqrtRhoInf(physical).';
            local = element.Volume.Weights ...
                * reshape(weighted, element.Volume.NPoints, []);
            determinants = mesh.computeElementJacobianDeterminants();
            expected = local*determinants(:);
            testCase.verifyEqual(weights.'*coefficients, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(space.eval(nodes, constant), ...
                ones(size(weighted)), 'AbsTol', 1e-12);
        end

        function vlasovPoissonMatchesPositiveLaplacian(testCase)
            mesh = mdemat.data.mesh.Grid({[0, 0.12, 0.6, 1]});
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 2);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            velocity = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1, 5));
            force = @(x) 0.1*cos(2*pi*x(1, :));
            penalty = {'all', 'right'; 'left', ''};
            config = struct('epsilon', 0.4, 'tau0', 2, ...
                'nh', 5, 'xPenaltyType', {penalty});
            config.ic = @(x, nh) [1+force(x); zeros(nh-1, size(x, 2))];
            kinetic = mdemat.application.plasma.HermiteState(space, velocity);
            scheme = mdemat.application.plasma.DghsScheme(config=config);
            scheme.initialize(kinetic);
            poisson = mdemat.solver.elliptic.LdgScheme(config= ...
                struct('force', force, 'xPenaltyType', {penalty}));
            potential = mdemat.common.fem.SpatialState(space);
            poisson.initialize(potential);
            poisson.step(potential);
            testCase.verifyEqual(kinetic.Dofs.P, potential.Dofs.U, 'AbsTol', 1e-11);
            testCase.verifyEqual(kinetic.Dofs.E, -potential.Dofs.Q, 'AbsTol', 1e-11);
            nodes = element.Volume.Nodes;
            local = element.Volume.Weights ...
                * reshape(space.eval(nodes, kinetic.Dofs.P), ...
                element.Volume.NPoints, []);
            determinants = mesh.computeElementJacobianDeterminants();
            testCase.verifyEqual(local*determinants(:), 0, 'AbsTol', 1e-12);
        end
    end

    methods (Static, Access = private)
        function [xSpace, vSpace, coefficients, multiplier] = fixture()
            mesh = mdemat.data.mesh.UniformGrid([0, 1], 2);
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 1);
            xSpace = mdemat.common.fem.FiniteElementSpace(mesh, element);
            hermite = mdemat.common.fem.L2OrthotopeElement.hermite(1, 5);
            vSpace = mdemat.common.fem.SpectralSpace(hermite);
            ng = xSpace.NGlobalDofs;
            coefficients = reshape(sin(1:5*ng), ng, 5);
            multiplier = diag(linspace(0.2, 0.8, ng)) ...
                + diag(ones(ng-1, 1)/7, -1);
        end

        function expected = blockSolve(initial, spatial, equilibrium, dt, method)
            nh = size(initial, 2);
            lower = diag(sqrt(1:nh-1), -1);
            generator = kron(lower, spatial);
            identity = eye(numel(initial));
            u = initial(:);
            e = equilibrium(:);
            switch method
                case 'be'
                    u = (identity-dt*generator) \ (u-dt*generator*e);
                case 'midpoint'
                    u = (identity-dt*generator/2) ...
                        \ ((identity+dt*generator/2)*u-dt*generator*e);
                case 'sdirk2'
                    gamma = 1-sqrt(2)/2;
                    stageMatrix = identity-gamma*dt*generator;
                    first = stageMatrix \ (u-gamma*dt*generator*e);
                    slope = generator*(first-e);
                    second = stageMatrix ...
                        \ (u+(1-gamma)*dt*slope-gamma*dt*generator*e);
                    u = u+dt*((1-gamma)*slope+gamma*generator*(second-e));
            end
            expected = reshape(u, size(initial));
        end
    end
end
