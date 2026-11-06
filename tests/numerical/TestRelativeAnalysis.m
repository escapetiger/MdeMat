classdef TestRelativeAnalysis < matlab.unittest.TestCase
    methods (Test)
        function explicitSameBasisPaddingMeasuresPhysicalArea(test)
            errors = zeros(2,3);
            componentErrors = zeros(2,3);
            scales = zeros(2,1);
            for level = 1:2
                [fine,coarse] = states(level);
                plan = mdemat.common.fem.nestedSamplingPlan(fine.XDisc,coarse.XDisc,3);
                values = fine.XDisc.eval(plan.FineNodes,fine.Dofs.U(plan.FineDofIndices,:)).';
                reference = coarse.XDisc.eval(plan.CoarseNodes,coarse.Dofs.U).';
                % These ordinary scalar components share an explicit identity;
                % component 3 is absent (zero) in the coarse representation.
                reference(3,:) = 0;
                report = mdemat.tools.analyze.weightedError(values,reference,plan.Weights);
                errors(level,:) = [norm(report.ComponentL1),norm(report.ComponentL2), ...
                    norm(report.ComponentLInf)];
                componentErrors(level,:) = report.ComponentL2.';
                scales(level) = coarse.XDisc.Mesh.computeMeasure();
            end
            magnitude = sqrt(0.4^2+0.1^2+0.3^2);
            test.verifyEqual(errors,[1;1/4]*[6,sqrt(6),1]*magnitude,AbsTol=2e-12);
            test.verifyEqual(componentErrors(:,[1,3]), ...
                sqrt(6)*[1;1/4]*[0.4,0.3],AbsTol=2e-12);
            rates = mdemat.tools.analyze.convergenceRates(scales,errors);
            test.verifyEqual(rates.LocalRates,2*ones(1,3),AbsTol=2e-12);
        end
    end
end

function [fine,coarse] = states(level)
element = mdemat.common.fem.BH1OrthotopeElement.nodal(2,3);
coarseSpace = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([0,2,-1,2],[1,2]*2^(level-1)),element);
fineSpace = coarseSpace.refine(1);
base = @(x) 1+x(1,:)+x(2,:).^2;
scale = 4^(1-level);
fine = struct('XDisc',fineSpace,'Dofs',struct('U', ...
    [fineSpace.fit(base),fineSpace.fit(0.2*scale),fineSpace.fit(0.3*scale)]));
coarse = struct('XDisc',coarseSpace,'Dofs',struct('U', ...
    [coarseSpace.fit(@(x) base(x)+0.4*scale),coarseSpace.fit(0.1*scale)]));
end
