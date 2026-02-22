classdef TestDynamicTimeline < matlab.unittest.TestCase
    
    properties (TestParameter)
        FinalTime = {1.0, 5.0, 10.0}
        NSteps = {1, 2, 5}
    end
    
    methods (Test)
        function testConstructor(testCase, FinalTime, NSteps)
            timeline = approx.mesh.DynamicTimeline(FinalTime, nSteps=NSteps);
            
            testCase.verifyEqual(timeline.Final, FinalTime);
            testCase.verifyEqual(timeline.Now, 0);
            testCase.verifyEqual(timeline.Count, 1);
            testCase.verifyEqual(timeline.Next, 0);
            testCase.verifyEqual(size(timeline.StepSizeQueue), [1, NSteps+1]);
            testCase.verifyEqual(size(timeline.OldStepSizeQueue), [1, NSteps+1]);
        end
        
        function testConstructorDefaults(testCase)
            timeline = approx.mesh.DynamicTimeline(5.0);
            
            testCase.verifyEqual(size(timeline.StepSizeQueue), [1, 2]);
        end
        
        function testConstructorValidation(testCase)                        
            testCase.verifyError(@() approx.mesh.DynamicTimeline(-1, nSteps=1), ...
                'MATLAB:validators:mustBeNonnegative');
        end
        
        function testSetTimeStep(testCase, FinalTime)
            timeline = approx.mesh.DynamicTimeline(FinalTime);
            
            h = 0.1;
            C = 0.5;
            p = 2;
            
            timeline.setTimeStep(h=h, C=C, p=p);
            
            expectedStepSize = C * h^p;
            testCase.verifyEqual(timeline.StepSizeQueue(1), expectedStepSize, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.Next, expectedStepSize, 'AbsTol', 1e-10);
            
            timeline.Now = FinalTime - 0.001;
            timeline.setTimeStep(h=h, C=C, p=p);
            
            testCase.verifyEqual(timeline.StepSizeQueue(1), FinalTime - timeline.Now, 'AbsTol', 1e-10);
        end
        
        function testSetTimeStepValidation(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            testCase.verifyError(@() timeline.setTimeStep(h=0.1, C=-0.5), ...
                'MATLAB:validators:mustBePositive');
        end
        
        function testAdvance(testCase, FinalTime)
            timeline = approx.mesh.DynamicTimeline(FinalTime);
            
            h = 0.1;
            C = 0.5;
            dt = C * h;
            timeline.setTimeStep(h=h, C=C);
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.Now, dt, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.Count, 2);
        end
        
        function testAdvanceMultiStep(testCase, FinalTime, NSteps)
            timeline = approx.mesh.DynamicTimeline(FinalTime, nSteps=NSteps);
            
            if NSteps >= 3
                timeline.StepSizeQueue = [0.1, 0.05, 0.025, 0.0125];
            else
                timeline.StepSizeQueue = 0.1 * ones(1, NSteps + 1);
            end
            
            initialH = timeline.StepSizeQueue;
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.StepSizeQueue, circshift(initialH, 1), 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.OldStepSizeQueue, initialH, 'AbsTol', 1e-10);
        end
        
        function testAdvanceToFinal(testCase, FinalTime)
            timeline = approx.mesh.DynamicTimeline(FinalTime);
            
            timeline.setTimeStep(h=FinalTime, C=1.0);
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.Now, FinalTime, 'AbsTol', 1e-10);
            testCase.verifyTrue(timeline.IsExhausted);
            
            testCase.verifyError(@() timeline.advance(), ...
                'approx:mesh:DynamicTimeline:TimeExceeded');
        end
        
        function testIsFinished(testCase, FinalTime)
            timeline = approx.mesh.DynamicTimeline(FinalTime);
            
            testCase.verifyFalse(timeline.IsExhausted);
            
            timeline.setTimeStep(h=FinalTime, C=1.0);
            
            testCase.verifyFalse(timeline.IsExhausted);
            
            timeline.advance();
            
            testCase.verifyTrue(timeline.IsExhausted);
        end
        
        function testHasStepSizeChanged(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            testCase.verifyTrue(timeline.HasStepSizeChanged);
            
            timeline.setTimeStep(h=0.1, C=1.0);
            
            timeline.advance();
            
            timeline.setTimeStep(h=0.05, C=1.0);
            
            testCase.verifyTrue(timeline.HasStepSizeChanged);
        end
        
        function testGetTimeStep(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            dt = 0.25;
            timeline.setTimeStep(h=dt, C=1.0);
            
            testCase.verifyEqual(timeline.StepSize, dt, 'AbsTol', 1e-10);
            
            dt2 = 0.5;
            timeline.setTimeStep(h=dt2, C=1.0);
            testCase.verifyEqual(timeline.StepSize, dt2, 'AbsTol', 1e-10);
        end
        
        function testMultiStepSequence(testCase, FinalTime)
            timeline = approx.mesh.DynamicTimeline(FinalTime, nSteps=2);
            
            testCase.verifyEqual(timeline.Now, 0);
            testCase.verifyEqual(timeline.Count, 1);
            
            timeline.setTimeStep(h=0.4, C=1.0, p=1);
            testCase.verifyEqual(timeline.StepSize, 0.16, 'AbsTol', 1e-10);
            
            timeline.advance();
            testCase.verifyEqual(timeline.Now, 0.16, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.Count, 2);
            
            timeline.setTimeStep(h=0.3, C=1.0, p=1);
            testCase.verifyTrue(timeline.HasStepSizeChanged);
            
            timeline.advance();
            testCase.verifyEqual(timeline.Now, 0.46, 'AbsTol', 1e-10);
            
            h = FinalTime-timeline.Now;
            timeline.setTimeStep(h=h, C=1.0);
            testCase.verifyEqual(timeline.StepSize, h, 'AbsTol', 1e-10);
            
            timeline.advance();
            testCase.verifyEqual(timeline.Now, FinalTime, 'AbsTol', 1e-10);
        end
        
        function testReset(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            timeline.setTimeStep(h=0.5, C=1.0);
            timeline.advance();
            
            timeline.reset();
            
            testCase.verifyEqual(timeline.Now, 0);
            testCase.verifyEqual(timeline.Count, 1);
        end
    end
end