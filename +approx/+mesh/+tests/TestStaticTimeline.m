classdef TestStaticTimeline < matlab.unittest.TestCase
    
    properties (TestParameter)
        TimeNodes = {[0, 0.1, 0.2, 0.5, 1.0], ...
                     [0, 0.25, 0.5, 0.75, 1.0], ...
                     [0, 1.0, 3.0, 6.0, 10.0]}
    end
    
    methods (Test)
        function testConstructor(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            testCase.verifyEqual(timeline.Final, TimeNodes(end));
            testCase.verifyEqual(timeline.Now, TimeNodes(1));
            testCase.verifyEqual(timeline.Count, 1);
            testCase.verifyEqual(timeline.Nodes, [TimeNodes, inf]);
        end
        
        function testConstructorValidation(testCase)
            testCase.verifyError(@() approx.mesh.StaticTimeline([0, 0.2, 0.1, 0.5]), ...
                'approx:mesh:StaticTimeline:InvalidInput');
        end
        
        function testAdvance(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            for i = 2:length(TimeNodes)
                timeline.advance();
                testCase.verifyEqual(timeline.Now, TimeNodes(i));
                testCase.verifyEqual(timeline.Count, i);
            end
        end
        
        function testIsFinished(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            testCase.verifyFalse(timeline.IsExhausted());
            
            for i = 2:length(TimeNodes)
                timeline.advance();
                testCase.verifyFalse(timeline.IsExhausted());
            end
            
            timeline.advance();
            testCase.verifyTrue(timeline.IsExhausted());
        end
        
        function testGetTimeStep(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            for i = 1:length(TimeNodes)-1
                expectedStepSize = TimeNodes(i+1) - TimeNodes(i);
                testCase.verifyEqual(timeline.StepSize, expectedStepSize, 'AbsTol', 1e-10);
                
                if i < length(TimeNodes)-1
                    timeline.advance();
                end
            end
            
            timeline.advance();
            
            testCase.verifyTrue(isinf(timeline.StepSize));
        end
        
        function testCountTotalTimeSteps(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            expectedCount = length(TimeNodes);
            testCase.verifyEqual(timeline.NTotalTimeSteps, expectedCount);
        end
        
        function testReset(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            for i = 1:2
                if i <= length(TimeNodes)-1
                    timeline.advance();
                end
            end
            
            timeline.reset();
            
            testCase.verifyEqual(timeline.Now, 0);
            testCase.verifyEqual(timeline.Count, 1);
        end
        
        function testAdvanceValidation(testCase, TimeNodes)
            timeline = approx.mesh.StaticTimeline(TimeNodes);
            
            for i = 1:length(TimeNodes)
                timeline.advance();
            end
            
            testCase.verifyError(@() timeline.advance(), ...
                'approx:mesh:StaticTimeline:TimeExceeded');
        end
    end
end