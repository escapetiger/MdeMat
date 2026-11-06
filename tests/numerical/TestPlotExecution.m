classdef TestPlotExecution < matlab.unittest.TestCase
    % Exercise sampling, controller lifetimes and real native graphics together.
    properties (TestParameter)
        family = {'elliptic','hyperbolic','parabolic','radiation','ap','dghs'}
    end
    methods (Test)
        function everyFamilyRendersItsActualFields(testCase,family)
            before = findall(groot,'Type','figure');
            session = mdemat.tools.visual.Visualizer(visible="off");
            testCase.addTeardown(@() session.finalize(force=true));
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([0,2*pi],3), ...
                mdemat.common.fem.BH1OrthotopeElement.nodal(1,2));
            config = struct('ic',@(x) ones(1,size(x,2)), ...
                'dt',0.001,'cfl',[],'tFinal',0.002);
            mode = '1d'; field = 'U';
            switch family
                case 'elliptic'
                    config = struct('bc',@(x) ones(1,size(x,2)), ...
                        'force',@(x) zeros(1,size(x,2)));
                    scheme = mdemat.solver.elliptic.PdgScheme(config=config);
                    state = mdemat.common.fem.SpatialState(space);
                case 'hyperbolic'
                    config.advection = 1;
                    scheme = mdemat.solver.hyperbolic.UwdgScheme(config=config);
                    state = mdemat.common.fem.SpatialState(space);
                case 'parabolic'
                    config.diffusion = 1;
                    scheme = mdemat.solver.parabolic.LdgScheme(config=config);
                    state = mdemat.common.fem.SpatialState(space);
                case 'radiation'
                    velocity = mdemat.common.fem.SumSpace( ...
                        mdemat.common.fem.SpectralSpace( ...
                        mdemat.common.fem.L2SphereElement.modal(1,1,reduction='topology')), ...
                        mdemat.common.fem.SpectralSpace( ...
                        mdemat.common.fem.L2SphereElement.nodal(1,2,reduction='topology')));
                    problem = mdemat.application.radiation.problem(space,velocity,1,0.002,0.001);
                    [result,scheme] = mdemat.application.radiation.solve(problem);
                    state = result.state; field = 'F';
                otherwise
                    velocity = mdemat.common.fem.SpectralSpace( ...
                        mdemat.common.fem.L2OrthotopeElement.hermite(1,3,T=1));
                    initial = @(x,nh) [ones(1,size(x,2));zeros(nh-1,size(x,2))];
                    problem = mdemat.application.plasma.problem(space,velocity,initial, ...
                        0.002,0.001,method=family,config=struct('vBBox',[-3,3]));
                    [result,scheme] = mdemat.application.plasma.solve(problem);
                    state = result.state; field = 'F'; mode = '2d';
            end
            scheme.setObserver('visualizer',session);
            scheme.setObservation('plot',struct('PointsPerDimension',2, ...
                'Components',struct(field,1),'Modes',{{mode}},'Times',[], ...
                'DatasetTypes',struct('NUMERIC','numeric')));
            state = scheme.run(state);
            testCase.verifyNotEmpty(fieldnames(session.FigureHandles));
            testCase.verifyTrue(session.IsEnabled);
            sampled = mdemat.internal.execution.sampleObservation(scheme,state,[]);
            values = sampled.NUMERIC.Data.(field).Values;
            testCase.verifyTrue(all(isfinite(values),'all'));
            if ~ismember(family,{'ap','dghs'})
                testCase.verifyEqual(values,ones(size(values)),'AbsTol',2e-10);
            else
                view = sampled.NUMERIC.Data.F;
                [~,v] = ndgrid(view.Axes{:});
                testCase.verifyEqual(values,exp(-v(:).'.^2/2)/sqrt(2*pi),'AbsTol',2e-10);
            end
            session.finalize(force=true);
            testCase.verifyEqual(findall(groot,'Type','figure'),before);
        end
    end
end
