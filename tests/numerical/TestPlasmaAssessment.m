classdef TestPlasmaAssessment < matlab.unittest.TestCase
    % Physical F is measured in dx dv on an explicitly truncated velocity domain.
    methods (Test)
        function hermiteAndNonconstantApScaleMatchPhysicalReferencePairs(testCase)
            space = spatialSpace([0,0.3,2]);
            time = 0.4;
            temperature = 0.7;
            maxwell = @(v) exp(-v.^2/(2*temperature))/sqrt(2*pi*temperature);
            for ap = [false,true]
                state = makeState(space,temperature,4,ap);
                if ap
                    base = 1;
                    scale = @(x) 2-0.25*x;
                    state.Coefs.sqrtRhoInf = space.fit(scale);
                else
                    base = 2;
                    scale = @(x) ones(size(x));
                end
                state.Dofs.D(:,1) = space.fit(@(x) base+x+time);
                scheme = makeScheme(ap,temperature,4,[-1,1]);
                references = { ...
                    @(x,v) scale(x).*(base+x+time).*maxwell(v), ...
                    @(x,v,t) scale(x).*(base+x+t).*maxwell(v), ...
                    @(x,v,varargin) scale(x).*(base+x+varargin{2}).*maxwell(v)};
                for index = 1:numel(references)
                    scheme.setObservation('analysis',struct('PointsPerDimension',3, ...
                        'Exacts',struct('F',references{index})));
                    sample = scheme.samplePhysicalAssessment(state,[],'F',time);
                    [X,V] = ndgrid(sample.SpatialCoordinates(1,:),sample.VelocityNodes);
                    expected = scale(X).*(base+X+time).*maxwell(V);
                    testCase.verifyEqual(sample.Values,expected(:).','AbsTol',5e-13);
                    testCase.verifyEqual(sample.Reference,expected(:).','AbsTol',5e-13);
                    testCase.verifyEqual(sum(sample.Weights),4,'AbsTol',5e-14);
                    testCase.verifyTrue(all(sample.Weights>0));
                    testCase.verifyEqual(sample.ComponentIds,1);
                    testCase.verifyEqual(sample.Comparison,'exact');
                    testCase.verifyTrue(contains(sample.ReferenceId,'truncated v='));
                    report = mdemat.tools.analyze.weightedError( ...
                        sample.Values,sample.Reference,sample.Weights);
                    testCase.verifyLessThan(report.L2,5e-13);
                end
            end
        end

        function sameSizeDifferentVelocityBasesAreComparedByPhysicalReconstruction(testCase)
            coarseSpace = spatialSpace([0,0.3,2]);
            fineSpace = coarseSpace.refine(1);
            fineTemperature = 0.8;
            coarseTemperature = 1.3;
            fine = makeState(fineSpace,fineTemperature,8,false);
            coarse = makeState(coarseSpace,coarseTemperature,8,false);
            fine.Dofs.D(:,1) = fineSpace.fit(1);
            coarse.Dofs.D(:,1) = coarseSpace.fit(1);
            scheme = makeScheme(false,fineTemperature,8,[-1,1]);
            scheme.setObservation('analysis',struct('PointsPerDimension',3));
            sample = scheme.samplePhysicalAssessment(fine,coarse,'F',[]);
            [~,V] = ndgrid(sample.SpatialCoordinates(1,:),sample.VelocityNodes);
            expectedFine = exp(-V.^2/(2*fineTemperature))/sqrt(2*pi*fineTemperature);
            expectedCoarse = exp(-V.^2/(2*coarseTemperature))/sqrt(2*pi*coarseTemperature);
            testCase.verifyEqual(sample.Values,expectedFine(:).','AbsTol',2e-13);
            testCase.verifyEqual(sample.Reference,expectedCoarse(:).','AbsTol',2e-13);
            report = mdemat.tools.analyze.weightedError(sample.Values,sample.Reference,sample.Weights);
            % Independent continuous Gaussian-product integral on x in [0,2], v in [-1,1].
            a = fineTemperature;
            b = coarseTemperature;
            square = @(T) erf(1/sqrt(T))/(2*sqrt(pi*T));
            cross = erf(sqrt((a+b)/(2*a*b)))/sqrt(2*pi*(a+b));
            expectedL2 = sqrt(2*(square(a)+square(b)-2*cross));
            testCase.verifyEqual(report.L2,expectedL2,'AbsTol',3e-13);
            testCase.verifyGreaterThan(report.L2,0.01);
            testCase.verifyEqual(sample.Comparison,'interlevel');
            testCase.verifyEqual(sample.VelocityBounds,[-1,1]);
        end

        function interlevelApUsesItsOwnNonconstantScaleOnReorderedSpatialPoints(testCase)
            coarseSpace = spatialSpace([0,0.3,2]);
            fineSpace = coarseSpace.refine(1);
            fine = makeState(fineSpace,0.9,4,true);
            coarse = makeState(coarseSpace,1.2,6,true);
            fine.Dofs.D(:,1) = fineSpace.fit(@(x) 1+x);
            coarse.Dofs.D(:,1) = coarseSpace.fit(@(x) 2-x/4);
            fine.Coefs.sqrtRhoInf = fineSpace.fit(@(x) 1+x/3);
            coarse.Coefs.sqrtRhoInf = coarseSpace.fit(@(x) 2-x/5);
            scheme = makeScheme(true,0.9,4,[-1,1]);
            scheme.setObservation('analysis',struct('PointsPerDimension',3));
            sample = scheme.samplePhysicalAssessment(fine,coarse,'F',[]);
            [X,V] = ndgrid(sample.SpatialCoordinates(1,:),sample.VelocityNodes);
            expectedFine = (1+X).*(1+X/3).*exp(-V.^2/(2*0.9))/sqrt(2*pi*0.9);
            expectedCoarse = (2-X/4).*(2-X/5).*exp(-V.^2/(2*1.2))/sqrt(2*pi*1.2);
            testCase.verifyEqual(sample.Values,expectedFine(:).','AbsTol',8e-13);
            testCase.verifyEqual(sample.Reference,expectedCoarse(:).','AbsTol',8e-13);
            testCase.verifyEqual(numel(sample.VelocityNodes),12);
        end

        function missingCoarseAndInvalidPhysicalReferencesFailExplicitly(testCase)
            space = spatialSpace([0,1]);
            state = makeState(space,1,4,false);
            state.Dofs.D(:,1) = space.fit(1);
            scheme = makeScheme(false,1,4,[-1,1]);
            scheme.setObservation('analysis',struct('PointsPerDimension',3));
            testCase.verifyEmpty(scheme.samplePhysicalAssessment(state,[],'F',0));
            testCase.verifyError(@() scheme.samplePhysicalAssessment(state,[],'D',0), ...
                'mdemat:plasma:assessmentSamples:Field');
            scheme.setObservation('analysis',struct('Exacts',struct('F',@(x,v) NaN)));
            testCase.verifyError(@() scheme.samplePhysicalAssessment(state,[],'F',0), ...
                'mdemat:plasma:assessmentSamples:NonFinite');
            scheme.setObservation('analysis',struct('Exacts',struct('F',@(x,v) zeros(2,numel(x)))));
            testCase.verifyError(@() scheme.samplePhysicalAssessment(state,[],'F',0), ...
                'mdemat:plasma:referenceDistribution:Shape');
            invalid = makeScheme(false,1,4,[-Inf,Inf]);
            invalid.setObservation('analysis',struct('PointsPerDimension',3,'Exacts',struct('F',0)));
            testCase.verifyError(@() invalid.samplePhysicalAssessment(state,[],'F',0), ...
                'mdemat:plasma:assessmentSamples:VelocityDomain');
        end
    end
end

function space = spatialSpace(vertices)
space = mdemat.common.fem.FiniteElementSpace(mdemat.data.mesh.Grid({vertices}), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(1,3));
end

function state = makeState(space,temperature,modes,ap)
velocity = mdemat.common.fem.SpectralSpace( ...
    mdemat.common.fem.L2OrthotopeElement.hermite(1,modes,T=temperature));
if ap
    state = mdemat.application.plasma.ApHermiteState(space,velocity);
    state.Coefs.sqrtRhoInf = space.fit(1);
else
    state = mdemat.application.plasma.HermiteState(space,velocity);
end
state.Coefs.T0 = temperature;
state.Dofs.D = zeros(space.NGlobalDofs,modes);
end

function scheme = makeScheme(ap,temperature,modes,bounds)
config = struct('dt',0.1,'cfl',[],'tFinal',0.2,'T0',temperature,'nh',modes,'vBBox',bounds);
if ap
    scheme = mdemat.application.plasma.ApDghsScheme(config=config);
else
    scheme = mdemat.application.plasma.DghsScheme(config=config);
end
end
