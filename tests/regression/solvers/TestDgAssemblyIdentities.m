classdef TestDgAssemblyIdentities < matlab.unittest.TestCase
    %TESTDGASSEMBLYIDENTITIES Verify weak derivatives without solving a PDE.
    methods (Test)
        function polynomialWeakDerivativesIncludeEachFaceExactlyOnce(testCase)
            for nd = 1:2
                for cells = 1:3
                    space = TestDgAssemblyIdentities.space(nd, cells);
                    assembly = mdemat.common.fem.AdjointAssembly(space, 'dirichlet', ...
                        {'boundary', 'right'; 'left', ''});
                    exact = @(x) sum(x.^2+0.3*x, 1)+0.7;
                    values = space.fit(exact);
                    matrix = assembly.assembleMatrix(eye(nd));
                    boundary = assembly.assembleVector(eye(nd), exact);
                    primal = zeros(space.NGlobalDofs, 1);
                    for d = 1:nd
                        derivative = space.fit(@(x) 2*x(d, :)+0.3);
                        testCase.verifyEqual(full(matrix{2, 1}{d}*values+boundary{2, 1}{d}), ...
                            full(derivative), 'AbsTol', 1e-11);
                        primal = primal-matrix{1, 1}{d}*values ...
                            -matrix{1, 2}{d}*derivative-boundary{1, 1}{d};
                    end
                    source = mdemat.common.fem.SourceAssembly(space);
                    testCase.verifyEqual(full(primal), full(source.assembleVector(-2*nd)), ...
                        'AbsTol', 1e-11);
                end
            end
        end

        function upwindPreservesPeriodicConstantsAndPhysicalMass(testCase)
            for nd = 1:2
                for cells = 1:3
                    space = TestDgAssemblyIdentities.space(nd, cells);
                    assembly = mdemat.common.fem.UpwindAssembly(space, 'periodic');
                    matrix = assembly.assembleMatrix((1:nd).');
                    constant = space.fit(1);
                    weights = space.integrationWeights();
                    testCase.verifyEqual(full(matrix*constant), zeros(size(constant)), 'AbsTol', 1e-11);
                    testCase.verifyEqual(full(weights.'*matrix), zeros(1, numel(constant)), 'AbsTol', 1e-11);
                end
            end
        end

        function constantWeightedAndPointwisePathsRespectAdjointSigns(testCase)
            for cells = 1:5
                space = TestDgAssemblyIdentities.space(1, cells);
                penalty = {'boundary', 'right'; 'left', ''};
                ordinary = mdemat.common.fem.AdjointAssembly(space, 'dirichlet', penalty);
                weighted = mdemat.common.fem.WeightedAdjointAssembly(space, 'dirichlet', penalty);
                expected = ordinary.assembleMatrix(1);
                byFunction = ordinary.assembleMatrix(@(x) ones(1, size(x, 2)));
                byWeight = weighted.assembleMatrix(1, @(x) ones(1, size(x, 2)), ...
                    dWeight=@(x) zeros(1, size(x, 2)));
                for row = 1:2
                    for column = 1:2
                        if row == 2 && column == 2, continue; end
                        testCase.verifyEqual(full(byFunction{row, column}{1}), full(expected{row, column}{1}), ...
                            'AbsTol', 1e-11);
                        % Weighted first row represents A*=-c*d/dx;
                        % ordinary first row represents +d/dx.
                        sign = 1;
                        if row == 1, sign = -1; end
                        testCase.verifyEqual(full(byWeight{row, column}{1}), sign*full(expected{row, column}{1}), ...
                            'AbsTol', 1e-11);
                    end
                end
            end
        end

        function boundaryPenaltyDampsTheWeakGradientNullMode(testCase)
            space = TestDgAssemblyIdentities.space(1, 2);
            element = space.Element;
            reference = element.Volume.Nodes;
            % The final-cell highest Legendre mode is orthogonal to every
            % test derivative. Homogeneous alternating traces make D21*z=0.
            local = element.Approximator.fit(element.Approximator.embed( ...
                (reference(:).^2-1/12)));
            mode = [zeros(space.NLocalDofs, 1); local];
            assembly = mdemat.common.fem.AdjointAssembly(space, 'dirichlet', ...
                {'boundary', 'right'; 'left', ''});
            dual = assembly.assembleMatrix(1);
            testCase.verifyEqual(full(dual{2, 1}{1}*mode), zeros(size(mode)), 'AbsTol', 1e-12);
            diffusion = mdemat.common.fem.EllipticAssembly(space, 'dirichlet', ...
                {'boundary', 'right'; 'left', ''});
            generator = diffusion.assembleMatrix(1);
            determinants = space.Mesh.computeElementJacobianDeterminants();
            mass = kron(diag(determinants), element.Approximator.Mass);
            testCase.verifyLessThan(mode.'*mass*generator*mode, -1e-6);
            symmetricPart = (mass*generator+generator.'*mass)/2;
            testCase.verifyLessThanOrEqual(max(eig(full(symmetricPart))), 1e-11);
        end

        function constantTensorCallbacksMatchNumericMatricesAndDirichletLoads(testCase)
            space = TestDgAssemblyIdentities.space(2, 2);
            tensor = [2, 0.3; 0.3, 1];
            callback = @(x) repmat(tensor(:), 1, size(x, 2));
            assembly = mdemat.common.fem.AdjointAssembly(space, 'dirichlet', ...
                {'boundary', 'right'; 'left', ''});
            expected = assembly.assembleMatrix(tensor);
            actual = assembly.assembleMatrix(callback);
            for row = 1:2
                for column = 1:2
                    if row == 2 && column == 2, continue; end
                    for d = 1:2
                        testCase.verifyEqual(full(actual{row, column}{d}), ...
                            full(expected{row, column}{d}), 'AbsTol', 1e-11);
                    end
                end
            end
            boundary = @(x) 1+0.2*sum(x, 1);
            expected = assembly.assembleVector(tensor, boundary);
            actual = assembly.assembleVector(callback, boundary);
            for row = 1:2
                for d = 1:2
                    testCase.verifyEqual(full(actual{row, 1}{d}), ...
                        full(expected{row, 1}{d}), 'AbsTol', 1e-11);
                end
            end
        end

        function upwindInflowLoadCompletesTheConservativeDerivative(testCase)
            for nd = 1:2
                for cells = 1:3
                    space = TestDgAssemblyIdentities.space(nd, cells);
                    field = @(x) 1+sum(x.^2+0.3*x, 1);
                    coefficients = space.fit(field);
                    assembly = mdemat.common.fem.UpwindAssembly(space, 'dirichlet');
                    source = mdemat.common.fem.SourceAssembly(space);
                    for orientation = [-1, 1]
                        velocity = orientation*(1:nd).';
                        derivative = @(x) sum(velocity.*(2*x+0.3), 1);
                        matrix = assembly.assembleMatrix(velocity);
                        lift = assembly.assembleVector(velocity, velocity, field);
                        testCase.verifyEqual(full(matrix*coefficients+lift), ...
                            full(source.assembleVector(derivative)), 'AbsTol', 1e-11);
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function space = space(nd, cells)
            mesh = mdemat.data.mesh.UniformGrid(repmat([0, 1], 1, nd), cells*ones(1, nd));
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(nd, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end
    end
end
