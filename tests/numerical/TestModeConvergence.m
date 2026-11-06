classdef TestModeConvergence < matlab.unittest.TestCase
    % Dynamic all-mode identities remain one curve as the velocity space grows.
    methods (Test)
        function textAndComplexSelectionsCannotBecomeNumericIndices(testCase)
            space = makeSpace(1);
            state = mdemat.common.fem.SpatialState(space);
            state.Dofs.U = zeros(space.NGlobalDofs,130);
            scheme = mdemat.solver.elliptic.PdgScheme();
            for selected = {'ab',"all"+["",""],[1+1i,2],[1,2;3,4]}
                scheme.setObservation('analysis',struct('PointsPerDimension',3, ...
                    'Components',struct('U',selected{1})));
                testCase.verifyError(@() mdemat.internal.execution.assessState( ...
                    scheme,state,[],[]),'mdemat:execution:assessState:Components');
            end
        end

        function nestedAllModesFollowAnIndependentQuadraticErrorLaw(testCase)
            spatial = makeSpace(1);
            velocity = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1,2,T=1));
            coarse = mdemat.application.plasma.HermiteState(spatial,velocity);
            coarse.Dofs.D = zeros(spatial.NGlobalDofs,2);
            coarse.Dofs.D(:,1) = spatial.fit(1);
            scheme = mdemat.application.plasma.DghsScheme();
            scheme.setObservation('analysis',struct('PointsPerDimension',3, ...
                'Components',struct('D','all'), ...
                'Reductions',struct('D',struct('x',{{'L2'}},'v','sequence'))));
            observations = [];
            for level = 1:3
                fine = coarse.refine(1);
                h = fine.XDisc.Mesh.computeMeasure();
                fine.Dofs.D = zeros(fine.NXDofs,fine.NVDofs);
                fine.Dofs.D(:,1) = fine.XDisc.fit(h^2);
                observation = mdemat.internal.execution.assessState(scheme,fine,coarse,[]);
                testCase.verifyNumElements(observation,1);
                testCase.verifyEqual(observation.ComponentIds,'all');
                testCase.verifyEqual(observation.ResolvedComponentIds,1:fine.NVDofs);
                testCase.verifyEqual(observation.Error,3*h^2,'AbsTol',1e-13);
                if isempty(observations)
                    observations = observation;
                else
                    observations(end+1) = observation; %#ok<AGROW>
                end
                coarse = fine;
            end
            summary = mdemat.tools.analyze.summarizeConvergence(observations);
            testCase.verifyNumElements(summary,1);
            testCase.verifyEqual(summary.LocalRates,[2;2],'AbsTol',2e-12);
            tables = mdemat.internal.execution.convergenceTables(summary);
            testCase.verifyEqual(tables.D_all_sequence.NComponents,[4;8;16]);
            testCase.verifyEqual(tables.D_all_sequence.L2_Error,3*[.5;.25;.125].^2,'AbsTol',1e-13);
            testCase.verifyEqual(tables.D_all_sequence.Properties.UserData.Identity.ComponentIds,'all');
        end

        function realApStudyRetainsBothObservedLevelsAndModeCounts(testCase)
            space = makeSpace(2);
            velocity = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1,3,T=1));
            config = struct('nh',3,'dt',0.001,'cfl',[],'tFinal',0.002);
            scheme = mdemat.application.plasma.ApDghsScheme(config=config);
            scheme.setObservation('analysis',struct('PointsPerDimension',3,'NLevels',2, ...
                'Components',struct('D','all'), ...
                'Reductions',struct('D',struct('x',{{'L2'}},'v','sequence'))));
            [state,tables] = scheme.converge(mdemat.application.plasma.ApHermiteState(space,velocity));
            testCase.verifyEqual(state.NVDofs,12);
            testCase.verifyEqual(state.XDisc.NMeshElements,8);
            testCase.verifyEqual(fieldnames(tables),{'D_all_sequence'});
            testCase.verifyEqual(tables.D_all_sequence.NComponents,[6;12]);
            testCase.verifyLessThan(max(tables.D_all_sequence.L2_Error),2e-11);
        end

        function longAndAmbiguousKeysKeepTheirFullDistinctIdentities(testCase)
            ids = 2*(1:(namelengthmax+5));
            a = observation('U1',23,'none',1);
            b = observation('U',123,'none',2);
            c = observation('large',ids,'sequence',3);
            summaries = mdemat.tools.analyze.summarizeConvergence([a,b,c]);
            tables = mdemat.internal.execution.convergenceTables(summaries);
            names = fieldnames(tables);
            testCase.verifyNumElements(names,3);
            testCase.verifyTrue(all(cellfun(@isvarname,names)));
            testCase.verifyTrue(all(cellfun(@numel,names)<=namelengthmax));
            recovered = cell(1,3);
            errors = zeros(1,3);
            for index = 1:3
                recovered{index} = tables.(names{index}).Properties.UserData.Identity;
                errors(index) = tables.(names{index}).L2_Error;
            end
            testCase.verifyEqual(errors,[1,2,3]);
            testCase.verifyEqual(recovered{1}.Field,'U1');
            testCase.verifyEqual(recovered{2}.Field,'U');
            testCase.verifyEqual(recovered{3}.ComponentIds,ids);
        end
    end
end

function space = makeSpace(cells)
space = mdemat.common.fem.FiniteElementSpace(mdemat.data.mesh.UniformGrid([0,1],cells), ...
    mdemat.common.fem.BH1OrthotopeElement.nodal(1,2));
end

function record = observation(field,ids,componentNorm,errorValue)
record = struct('Scale',1,'Field',field,'ComponentIds',ids,'SpatialNorm','L2', ...
    'ComponentNorm',componentNorm,'Comparison','exact','Reference','independent', ...
    'Error',errorValue);
end
