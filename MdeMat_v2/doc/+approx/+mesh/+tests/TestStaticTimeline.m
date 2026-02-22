classdef TestStaticTimeline < matlab.unittest.TestCase
    
    properties (TestParameter)
        timeNodes = {[0, 0.1, 0.2, 0.5, 1.0], ...
                     [0, 0.25, 0.5, 0.75, 1.0], ...
                     [0, 1.0, 3.0, 6.0, 10.0]}
    end
    
    methods (Test)
        function testConstructor(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            testCase.verifyEqual(timeline.final, timeNodes(end));
            testCase.verifyEqual(timeline.now, timeNodes(1));
            testCase.verifyEqual(timeline.count, 1);
            testCase.verifyEqual(timeline.nodes, [timeNodes, inf]);
        end
        
        function testConstructorValidation(testCase)
            testCase.verifyError(@() approx.mesh.StaticTimeline([0, 0.2, 0.1, 0.5]), ...
                'approx:mesh:StaticTimeline:InvalidInput');
        end
        
        function testAdvance(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            for i = 2:length(timeNodes)
                timeline.advance();
                testCase.verifyEqual(timeline.now, timeNodes(i));
                testCase.verifyEqual(timeline.count, i);
            end
        end
        
        function testIsFinished(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            testCase.verifyFalse(timeline.isFinished());
            
            for i = 2:length(timeNodes)
                timeline.advance();
                testCase.verifyFalse(timeline.isFinished());
            end
            
            timeline.advance();
            testCase.verifyTrue(timeline.isFinished());
        end
        
        function testGetTimeStep(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            for i = 1:length(timeNodes)-1
                expectedDt = timeNodes(i+1) - timeNodes(i);
                testCase.verifyEqual(timeline.dt, expectedDt, 'AbsTol', 1e-10);
                
                if i < length(timeNodes)-1
                    timeline.advance();
                end
            end
            
            timeline.advance();
            
            testCase.verifyTrue(isinf(timeline.dt));
        end
        
        function testCountTotalTimeSteps(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            expectedCount = length(timeNodes);
            testCase.verifyEqual(timeline.nTotalTimeSteps, expectedCount);
        end
        
        function testReset(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            for i = 1:2
                if i <= length(timeNodes)-1
                    timeline.advance();
                end
            end
            
            timeline.reset();
            
            testCase.verifyEqual(timeline.now, 0);
            testCase.verifyEqual(timeline.count, 1);
        end
        
        function testAdvanceValidation(testCase, timeNodes)
            timeline = approx.mesh.StaticTimeline(timeNodes);
            
            for i = 1:length(timeNodes)
                timeline.advance();
            end
            
            testCase.verifyError(@() timeline.advance(), ...
                'approx:mesh:StaticTimeline:TimeExceeded');
        end
    end
end