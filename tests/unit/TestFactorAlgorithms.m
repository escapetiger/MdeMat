classdef TestFactorAlgorithms < matlab.unittest.TestCase
    methods (Test)
        function physicalQrAndCompressionBound(test)
            Mx = diag([1, 2, 3, 4]); Mv = diag([2, 3, 4]);
            U = diag(1 ./ sqrt(diag(Mx))); U = U(:, 1:3);
            V = diag(1 ./ sqrt(diag(Mv)));
            S = diag([3, 1, 1e-5]);
            [F, report] = mdemat.common.lowrank.compress(U, S, V, ...
                LeftMetric=Mx, RightMetric=Mv, RelativeTolerance=1e-4);
            test.verifyEqual(F.Rank, 2);
            test.verifyTrue(report.ToleranceMet);
            difference = U*S*V' - F.materialize();
            actualError = norm(chol(Mx)*difference*chol(Mv)', 'fro');
            test.verifyEqual(actualError, report.DiscardedNorm, AbsTol=2e-14);
            test.verifyEqual(F.U'*Mx*F.U, eye(2), AbsTol=2e-14);
            test.verifyEqual(F.V'*Mv*F.V, eye(2), AbsTol=2e-14);
            test.verifyEqual(F.weightedNorm(), sqrt(10), AbsTol=2e-14);
            sparseFactors = mdemat.common.lowrank.compress(speye(4), ...
                sparse(diag([3,2,1,0])), speye(4), RelativeTolerance=1e-12);
            test.verifyEqual(full(sparseFactors.materialize()), diag([3,2,1,0]), AbsTol=2e-14);
            test.verifyEqual(sparseFactors.Rank,3);
        end

        function rankLimitIsReportedHonestly(test)
            [~, report] = mdemat.common.lowrank.compress(eye(3), diag([3, 2, 1]), ...
                eye(3), RelativeTolerance=1e-12, MaxRank=1);
            test.verifyFalse(report.ToleranceMet);
            test.verifyTrue(report.RankLimitReached);
            test.verifyEqual(report.DiscardedNorm, sqrt(5), AbsTol=1e-14);
            test.verifyError(@() mdemat.common.lowrank.compress(eye(3), eye(3), ...
                eye(3), MaxRank=1.5), 'mdemat:lowrank:RankBounds');
        end

        function complexGaugeAndSampling(test)
            U = [1, 0; 0, 1i; 1, -1] / 2;
            V = [1, 1i; 2, -1; 0, 1] / 3;
            S = [2, 1i; 0.3, -1];
            Q = [1, 1i; 1i, 1] / sqrt(2);
            F = mdemat.data.tensor.MatrixFactors(U, S, V);
            G = F.withFactors(U*Q, Q'*S*Q, V*Q);
            test.verifyEqual(G.materialize(), F.materialize(), AbsTol=3e-15);
            test.verifyEqual(F.weightedNorm(), norm(F.materialize(), 'fro'), AbsTol=2e-15);
            full = F.materialize();
            test.verifyEqual(F.sample([1, 3], [2, 3]), full([1, 3], [2, 3]));
        end

        function separatedActionMatchesDenseReference(test)
            F = mdemat.data.tensor.MatrixFactors([1; 2; 3], 2, [1i; 2]);
            A = {diag([1, 2, 3]), [0, 1, 0; -1, 0, 0; 0, 0, 0]};
            B = {eye(2), [0, 1i; 2, 1]};
            expected = A{1}*F.materialize()*B{1}' + A{2}*F.materialize()*B{2}';
            actual = mdemat.common.operator.applySeparated(A, B, F);
            test.verifyEqual(actual.materialize(), expected, AbsTol=5e-13);
            zero = mdemat.common.operator.applySeparated({}, {}, F);
            test.verifyEqual(zero.Shape, F.Shape);
            test.verifyEqual(zero.Rank, 0);
            test.verifyEqual(zero.weightedNorm(), 0);
        end
    end
end
