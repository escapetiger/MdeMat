classdef TestSpatialAssessment < matlab.unittest.TestCase
    methods (Test)
        function physicalNormsIncludeNonuniformJacobianAndAllComponents(test)
            mesh = mdemat.data.mesh.Grid({[-1,-0.7,0.2,2]});
            space = mdemat.common.fem.FiniteElementSpace(mesh, ...
                mdemat.common.fem.BH1OrthotopeElement.modal(1,3));
            plan = mdemat.common.fem.samplingPlan(space,8);
            reference = [plan.Coordinates;plan.Coordinates.^2];
            measured = mdemat.tools.analyze.weightedError( ...
                reference+[1;-2],reference,plan.Weights);
            test.verifyEqual([measured.ComponentL1,measured.ComponentL2,measured.ComponentLInf].', ...
                [3,6;sqrt(3),2*sqrt(3);1,2],AbsTol=2e-13);
            test.verifyFalse(any(structfun(@isobject,plan)));
        end

        function nestedComparisonAlignsActualFineAndCoarseCoordinates(test)
            for dimensions = 1:3
                coarse = makeSpace(dimensions,2);
                for ratio = [1,2,3]
                    fine = makeSpace(dimensions,2*ratio);
                    plan = mdemat.common.fem.nestedSamplingPlan(fine,coarse,3);
                    polynomial = @(x) [sum(x.^2,1);sum(x,1)];
                    C = coarse.fit(polynomial);
                    F = fine.fit(@(x) polynomial(x)+[0.5;-0.75]);
                    fineValues = fine.eval(plan.FineNodes,F(plan.FineDofIndices,:)).';
                    coarseValues = coarse.eval(plan.CoarseNodes,C).';
                    actualFine = fine.Mesh.collocate(plan.FineNodes,plan.FineElementIndices);
                    test.verifyEqual(actualFine,plan.Coordinates,AbsTol=2e-14);
                    measured = mdemat.tools.analyze.weightedError(fineValues,coarseValues,plan.Weights);
                    test.verifyEqual([measured.ComponentL1,measured.ComponentL2,measured.ComponentLInf].', ...
                        repmat([0.5,0.75],3,1),AbsTol=5e-12);
                end
            end
        end

        function explicitObservationsRetainComponentIdentityAndRates(test)
            observations = struct([]);
            for level = 1:3
                count = 2^level;
                space = makeSpace(1,count);
                plan = mdemat.common.fem.samplingPlan(space,8);
                coefficients = space.fit(@(x) [ones(1,size(x,2));2*ones(1,size(x,2))]/count^2);
                values = mdemat.common.fem.sampleField(space,coefficients,plan);
                report = mdemat.tools.analyze.weightedError(values,zeros(size(values)),plan.Weights);
                for component = 1:2
                    item = observation(1/count,component,report.ComponentL2(component));
                    if isempty(observations), observations = item;
                    else, observations(end+1) = item; end
                end
            end
            summaries = mdemat.tools.analyze.summarizeConvergence(observations);
            test.verifyEqual(summaries(1).Errors,(1./[2;4;8]).^2,AbsTol=1e-14);
            test.verifyEqual(summaries(2).Errors,2*summaries(1).Errors,AbsTol=1e-14);
            test.verifyEqual(summaries(1).LocalRates,[2;2],AbsTol=2e-13);
            test.verifyEqual(summaries(2).Identity.ComponentIds,2);
        end

        function sequenceAndGlobalNormsAreDifferentMathematicalQuestions(test)
            report = mdemat.tools.analyze.weightedError([1,0;0,1],zeros(2),[0.5,0.5]);
            test.verifyEqual(report.L1,1);
            test.verifyEqual(norm(report.ComponentL1),1/sqrt(2),AbsTol=4*eps);
            test.verifyEqual(report.LInf,1);
            test.verifyEqual(norm(report.ComponentLInf),sqrt(2),AbsTol=4*eps);
            item = observation(0.5,[2,4],norm(report.ComponentL2));
            item.ComponentNorm = 'sequence';
            summary = mdemat.tools.analyze.summarizeConvergence(item);
            test.verifyEqual(summary.Identity.ComponentIds,[2,4]);
            item.SpatialNorm = 'TYPO';
            test.verifyError(@() mdemat.tools.analyze.summarizeConvergence(item), ...
                'mdemat:analysis:InvalidObservation');
        end

        function mismatchedMeshesAndAxisCountsFailClearly(test)
            coarse = makeSpace(2,2);
            test.verifyError(@() mdemat.common.fem.samplingPlan(coarse,[2,3,4]), ...
                'mdemat:analysis:InvalidSampling');
            fine = makeSpace(2,3);
            test.verifyError(@() mdemat.common.fem.nestedSamplingPlan(fine,coarse,2), ...
                'mdemat:analysis:NonNestedMesh');
        end

        function rowAndColumnBoundsConstructTheSameNestedPlan(test)
            coarse = makeSpace(2,2);
            fine = makeSpace(2,4);
            bounds = [-0.5,0.5,-0.5,0.5];
            fine.Element.Geometry = mdemat.data.geometry.Orthotope(bounds);
            coarse.Element.Geometry = mdemat.data.geometry.Orthotope(bounds);
            rowPlan = mdemat.common.fem.nestedSamplingPlan(fine,coarse,[2,3]);
            fine.Element.Geometry = mdemat.data.geometry.Orthotope(bounds.');
            coarse.Element.Geometry = mdemat.data.geometry.Orthotope(bounds.');
            columnPlan = mdemat.common.fem.nestedSamplingPlan(fine,coarse,[2,3]);
            test.verifyEqual(rowPlan,columnPlan);
        end
    end
end

function space = makeSpace(dimensions,count)
space = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid(repmat([0,1],1,dimensions),count*ones(1,dimensions)), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(dimensions,3));
end

function item = observation(scale,components,value)
item = struct('Scale',scale,'Field','U','ComponentIds',components,'SpatialNorm','L2', ...
    'ComponentNorm','none','Comparison','exact','Reference','analytic zero','Error',value);
end
