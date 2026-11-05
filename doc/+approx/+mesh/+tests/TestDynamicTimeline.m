classdef TestDynamicTimeline < matlab.unittest.TestCase
    
    properties (TestParameter)
        finalTime = {1.0, 5.0, 10.0}
        nSteps = {1, 2, 5}
    end
    
    methods (Test)
        function testConstructor(testCase, finalTime, nSteps)
            timeline = approx.mesh.DynamicTimeline(finalTime, nSteps);
            
            testCase.verifyEqual(timeline.final, finalTime);
            testCase.verifyEqual(timeline.now, 0);
            testCase.verifyEqual(timeline.count, 1);
            testCase.verifyEqual(timeline.next, 0);
            testCase.verifyEqual(size(timeline.h), [1, nSteps+1]);
            testCase.verifyEqual(size(timeline.h0), [1, nSteps+1]);
        end
        
        function testConstructorDefaults(testCase)
            timeline = approx.mesh.DynamicTimeline(5.0);
            
            testCase.verifyEqual(size(timeline.h), [1, 2]);
        end
        
        function testConstructorValidation(testCase)                        
            testCase.verifyError(@() approx.mesh.DynamicTimeline(-1, 1), ...
                'approx:mesh:Timeline:InvalidInput');
        end
        
        function testSetTimeStep(testCase, finalTime)
            timeline = approx.mesh.DynamicTimeline(finalTime);
            
            h = 0.1;
            C = 0.5;
            p = 2;
            
            timeline.setTimeStep(h, C, p);
            
            expectedDt = C * h^p;
            testCase.verifyEqual(timeline.h(1), expectedDt, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.next, expectedDt, 'AbsTol', 1e-10);
            
            timeline.now = finalTime - 0.001;
            timeline.setTimeStep(h, C, p);
            
            testCase.verifyEqual(timeline.h(1), finalTime - timeline.now, 'AbsTol', 1e-10);
        end
        
        function testSetTimeStepValidation(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            testCase.verifyError(@() timeline.setTimeStep(0.1, -0.5), ...
                'approx:mesh:DynamicTimeline:InvalidInput');
        end
        
        function testAdvance(testCase, finalTime)
            timeline = approx.mesh.DynamicTimeline(finalTime);
            
            h = 0.1;
            C = 0.5;
            dt = C * h;
            timeline.setTimeStep(h, C);
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.now, dt, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.count, 2);
        end
        
        function testAdvanceMultiStep(testCase, finalTime, nSteps)
            timeline = approx.mesh.DynamicTimeline(finalTime, nSteps);
            
            if nSteps >= 3
                timeline.h = [0.1, 0.05, 0.025, 0.0125];
            else
                timeline.h = 0.1 * ones(1, nSteps + 1);
            end
            
            initialH = timeline.h;
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.h, circshift(initialH, 1), 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.h0, initialH, 'AbsTol', 1e-10);
        end
        
        function testAdvanceToFinal(testCase, finalTime)
            timeline = approx.mesh.DynamicTimeline(finalTime);
            
            timeline.setTimeStep(finalTime, 1.0);
            
            timeline.advance();
            
            testCase.verifyEqual(timeline.now, finalTime, 'AbsTol', 1e-10);
            testCase.verifyTrue(timeline.isFinished());
            
            testCase.verifyError(@() timeline.advance(), ...
                'approx:mesh:DynamicTimeline:TimeExceeded');
        end
        
        function testIsFinished(testCase, finalTime)
            timeline = approx.mesh.DynamicTimeline(finalTime);
            
            testCase.verifyFalse(timeline.isFinished());
            
            timeline.setTimeStep(finalTime, 1.0);
            
            testCase.verifyFalse(timeline.isFinished());
            
            timeline.advance();
            
            testCase.verifyTrue(timeline.isFinished());
        end
        
        function testHasStepSizeChanged(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            testCase.verifyTrue(timeline.hasStepSizeChanged());
            
            timeline.setTimeStep(0.1, 1.0);
            
            timeline.advance();
            
            originalStep = timeline.h(1);
            timeline.setTimeStep(0.05, 1.0);
            
            testCase.verifyTrue(timeline.hasStepSizeChanged());
        end
        
        function testGetTimeStep(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            dt = 0.25;
            timeline.setTimeStep(dt, 1.0);
            
            testCase.verifyEqual(timeline.dt, dt, 'AbsTol', 1e-10);
            
            dt2 = 0.5;
            timeline.setTimeStep(dt2, 1.0);
            testCase.verifyEqual(timeline.dt, dt2, 'AbsTol', 1e-10);
        end
        
        function testMultiStepSequence(testCase, finalTime)
            timeline = approx.mesh.DynamicTimeline(finalTime, 2);
            
            testCase.verifyEqual(timeline.now, 0);
            testCase.verifyEqual(timeline.count, 1);
            
            timeline.setTimeStep(0.4, 1.0, 1);
            testCase.verifyEqual(timeline.dt, 0.4, 'AbsTol', 1e-10);
            
            timeline.advance();
            testCase.verifyEqual(timeline.now, 0.4, 'AbsTol', 1e-10);
            testCase.verifyEqual(timeline.count, 2);
            
            timeline.setTimeStep(0.3, 1.0, 1);
            testCase.verifyTrue(timeline.hasStepSizeChanged);
            
            timeline.advance();
            testCase.verifyEqual(timeline.now, 0.7, 'AbsTol', 1e-10);
            
            timeline.setTimeStep(0.3, 1.0);
            expectedDt = min(0.3, finalTime - timeline.now);
            testCase.verifyEqual(timeline.dt, expectedDt, 'AbsTol', 1e-10);
            
            timeline.advance();
            expectedFinal = min(1.0, finalTime);
            testCase.verifyEqual(timeline.now, expectedFinal, 'AbsTol', 1e-10);
        end
        
        function testReset(testCase)
            timeline = approx.mesh.DynamicTimeline(10.0);
            
            timeline.setTimeStep(0.5, 1.0);
            timeline.advance();
            
            timeline.reset();
            
            testCase.verifyEqual(timeline.now, 0);
            testCase.verifyEqual(timeline.count, 1);
        end
    end
end