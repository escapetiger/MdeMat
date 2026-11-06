classdef TestLowRankDynamics < matlab.unittest.TestCase
    methods (Test)
        function exactRankEvolutionAndSmallSingularValue(test)
            A = [0, 1, 0; -1, -0.2, 0; 0, 0, -0.4];
            B = [0, 0.3, 0, 0; -0.3, 0, 0, 0; 0, 0, -0.1, 0; 0, 0, 0, -0.2];
            F = mdemat.data.tensor.MatrixFactors(eye(3, 2), diag([1, 1e-12]), eye(4, 2));
            rhs = @(~, field) mdemat.common.operator.applySeparated({A, eye(3)}, {eye(4), B}, field);
            [actual, report] = mdemat.common.ode.evolveFixedRank(F, rhs, [0, 0.2], 0.01);
            expected = expm(0.2*A)*F.materialize()*expm(0.2*B)';
            test.verifyLessThan(norm(actual.materialize()-expected, 'fro'), 2e-10);
            test.verifyEqual(actual.Rank, 2);
            test.verifyEqual(report.SplittingOrder, 1);
        end

        function projectedFlowHasFirstOrderTimeRefinement(test)
            F0 = [1, 0; 0.2, 0; 0, 0];
            forcing = [0.1, 0.5; 0, -0.4; 0.6, 0.2];
            F = mdemat.data.tensor.MatrixFactors(F0(:, 1), 1, [1; 0]);
            G = mdemat.data.tensor.MatrixFactors(forcing, eye(2), eye(2));
            rhs = @(~, ~) G;
            [~, reference] = ode45(@(~, y) projected(y, forcing), [0, 0.3], F0(:), ...
                odeset('RelTol', 2e-11, 'AbsTol', 1e-13));
            expected = reshape(reference(end, :), 3, 2);
            errors = zeros(1, 3);
            for j = 1:3
                actual = mdemat.common.ode.evolveFixedRank(F, rhs, [0, 0.3], 0.03/2^(j-1));
                errors(j) = norm(actual.materialize()-expected, 'fro');
            end
            rates = log2(errors(1:end-1)./errors(2:end));
            test.verifyGreaterThan(rates, 0.9);
            test.verifyLessThan(rates, 1.15);
        end

        function adaptiveGrowthDropAndPortableRestart(test)
            F = mdemat.data.tensor.MatrixFactors([1; 0; 0], 1, [1; 0]);
            G = mdemat.data.tensor.MatrixFactors([0; 1; 0], 1, [0; 1]);
            [grown, growth] = mdemat.common.ode.evolveAdaptive(F, @(~, ~) G, [0, 0.2], 0.1, ...
                RelativeTolerance=1e-12);
            test.verifyEqual(grown.Rank, 2);
            test.verifyEqual(grown.materialize(), F.materialize()+0.2*G.materialize(), AbsTol=2e-14);
            snapshot = grown.snapshot();
            restored = mdemat.data.tensor.MatrixFactors.fromSnapshot(snapshot);
            [continued, ~] = mdemat.common.ode.evolveAdaptive(restored, @(~, ~) G, [0.2, 0.4], 0.1);
            direct = mdemat.common.ode.evolveAdaptive(F, @(~, ~) G, [0, 0.4], 0.1);
            test.verifyEqual(continued.materialize(), direct.materialize(), AbsTol=2e-14);
            [dropped, ~] = mdemat.common.ode.evolveAdaptive(grown, @(~, ~) G, [0.2, 0.3], 0.1, ...
                RelativeTolerance=0.4);
            test.verifyEqual(dropped.Rank, 1);
            test.verifyGreaterThan(max(growth.Ranks), min(growth.Ranks));
            snapshot.Version = 99;
            test.verifyError(@() mdemat.data.tensor.MatrixFactors.fromSnapshot(snapshot), ...
                'mdemat:tensor:SnapshotSchema');
        end

        function slabTransportMatchesDenseAndPreservesMass(test)
            nx = 12; v = [-1, -0.4, 0.3, 1]; w = [0.2; 0.4; 0.4; 0.2];
            h = 1/nx; x = ((1:nx)'-0.5)*h;
            F = mdemat.data.tensor.MatrixFactors([ones(nx, 1), sin(2*pi*x)], ...
                diag([1, 0.1]), [ones(4, 1), v'], ...
                LeftMetric=h*speye(nx), RightMetric=diag(w));
            duration = 0.04; scattering = 0.7;
            i = (1:nx)';
            Db = (eye(nx)-full(sparse(i, mod(i-2, nx)+1, 1, nx, nx)))/h;
            Df = (full(sparse(i, mod(i, nx)+1, 1, nx, nx))-eye(nx))/h;
            C = w*ones(1, 4)/sum(w)-eye(4);
            operator = -kron(diag(max(v, 0)), Db)-kron(diag(min(v, 0)), Df) ...
                +scattering*kron(C', eye(nx));
            reference = reshape(expm(duration*operator)*reshape(F.materialize(), [], 1), nx, 4);
            errors = zeros(1, 3);
            for j = 1:3
                [solution, report] = mdemat.application.radiation.solveLowRankTransport( ...
                    F, v, w, [0, duration], 0.004/2^(j-1), ...
                    Scattering=scattering, RelativeTolerance=1e-12);
                errors(j) = norm(solution.Field.materialize()-reference, 'fro');
                test.verifyEqual(report.FinalMass, report.InitialMass, AbsTol=5e-13);
            end
            test.verifyGreaterThan(log2(errors(1:end-1)./errors(2:end)), 0.95);
            [solution, report] = mdemat.application.radiation.solveLowRankTransport( ...
                F, v, w, [0, duration], 0.001, Scattering=scattering, ...
                MaxRank=2, RelativeTolerance=1e-4);
            test.verifyEqual(report.FinalMass, report.InitialMass, AbsTol=5e-13);
            test.verifyLessThanOrEqual(solution.Field.Rank, 2);
            test.verifyFalse(report.PositivityGuaranteedAfterCompression);
        end

        function kineticCheckpointRestartsTheSameAdaptiveTrajectory(test)
            nx = 10; nv = 4; v = [-1,-0.3,0.3,1]; w = ones(nv,1)/nv;
            h = 1/nx; x = ((1:nx)'-0.5)*h;
            F = mdemat.data.tensor.MatrixFactors([ones(nx,1),cos(2*pi*x)], ...
                diag([1,0.1]),[ones(nv,1),v'], ...
                LeftMetric=h*speye(nx),RightMetric=diag(w));
            fullRun = mdemat.application.radiation.solveLowRankTransport( ...
                F,v,w,[0,0.04],0.002,Scattering=0.5,MaxRank=2);
            halfRun = mdemat.application.radiation.solveLowRankTransport( ...
                F,v,w,[0,0.02],0.002,Scattering=0.5,MaxRank=2);
            file = string(tempname)+".mat";
            fileGuard = onCleanup(@() delete(file));
            record = struct('Field',halfRun.Field.snapshot(),'Time',halfRun.Time, ...
                'Step',0.002,'Velocities',v,'Weights',w,'Scattering',0.5);
            mdemat.tools.io.writeRecord(file,record);
            restored = mdemat.tools.io.readRecord(file);
            restart = mdemat.data.tensor.MatrixFactors.fromSnapshot(restored.Field);
            continued = mdemat.application.radiation.solveLowRankTransport( ...
                restart,restored.Velocities,restored.Weights,[restored.Time,0.04],restored.Step, ...
                Scattering=restored.Scattering,MaxRank=2);
            test.verifyEqual(continued.Field.materialize(),fullRun.Field.materialize(),AbsTol=5e-13);
            test.verifyGreaterThan(min(continued.Field.materialize(),[],'all'),0);
        end
    end
end

function value = projected(y, forcing)
Y = reshape(y, 3, 2);
[U, ~, V] = svd(Y, 'econ');
U = U(:, 1); V = V(:, 1);
value = U*(U'*forcing)+(forcing*V)*V'-U*(U'*forcing*V)*V';
value = value(:);
end
