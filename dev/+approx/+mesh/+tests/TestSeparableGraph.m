classdef TestSeparableGraph < matlab.unittest.TestCase
    
    methods (Test)
        function testConstructor(testCase)
            V1D = {[0; 1; 2; 3]};
            E1D = {[1, 2; 2, 3; 3, 4]};
            graph1D = approx.mesh.SeparableGraph(V1D, E1D);
            testCase.verifyEqual(graph1D.Vertices, V1D);
            testCase.verifyEqual(graph1D.Edges, E1D);
            
            V2D = {[0; 1; 2], [0; 1; 2]};
            E2D = {[1, 2; 2, 3], [1, 2; 2, 3]};
            graph2D = approx.mesh.SeparableGraph(V2D, E2D);
            testCase.verifyEqual(graph2D.Vertices, V2D);
            testCase.verifyEqual(graph2D.Edges, E2D);
            
            V3D = {[0; 1], [0; 1], [0; 1]};
            E3D = {[1, 2], [1, 2], [1, 2]};
            graph3D = approx.mesh.SeparableGraph(V3D, E3D);
            testCase.verifyEqual(graph3D.Vertices, V3D);
            testCase.verifyEqual(graph3D.Edges, E3D);
        end
        
        function testConstructorValidation(testCase)
            V = {[0; 1; 2]};
            E = {[1, 2; 2, 3]};
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(), ...
                'MATLAB:minrhs');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V), ...
                'MATLAB:minrhs');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph('not cell', E), ...
                'MATLAB:validation:UnableToConvert');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V, 'not cell'), ...
                'MATLAB:validation:UnableToConvert');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V, {[1, 2], [1, 2]}), ...
                'approx:mesh:SeparableGraph:DimensionMismatch');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph({[1, 2; 3, 4]}, E), ...
                'approx:mesh:SeparableGraph:InvalidVertices');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V, {[1, 2, 3]}), ...
                'approx:mesh:SeparableGraph:InvalidEdges');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V, {[1, 4]}), ...
                'approx:mesh:SeparableGraph:InvalidEdgeIndices');
        end
        
        function testNDims(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2]}, {[1, 2; 2, 3]});
            testCase.verifyEqual(graph1D.NDims, 1);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1], [0; 1]}, {[1, 2], [1, 2]});
            testCase.verifyEqual(graph2D.NDims, 2);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.NDims, 3);
        end
        
        function testNVertices(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.NVertices, 4);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.NVertices, [3, 3]);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.NVertices, [2, 2, 2]);
        end
        
        function testNTotalVertices(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.NTotalVertices, 4);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.NTotalVertices, 9);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.NTotalVertices, 8);
        end
        
        function testNEdges(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.NEdges, 3);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.NEdges, [2, 2]);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.NEdges, [1, 1, 1]);
        end
        
        function testNTotalEdges(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.NTotalEdges, 3);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.NTotalEdges, 12);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.NTotalEdges, 12);
        end
        
        function testFull1D(testCase)
            V = {[0; 1; 2; 3]};
            E = {[1, 2; 2, 3; 3, 4]};
            graph = approx.mesh.SeparableGraph(V, E);
            fullGraph = graph.full();
            
            testCase.verifyClass(fullGraph, 'approx.mesh.Graph');
            testCase.verifyEqual(fullGraph.NVertices, 4);
            testCase.verifyEqual(fullGraph.NEdges, 3);
            testCase.verifyEqual(fullGraph.Vertices, V{1});
            testCase.verifyEqual(fullGraph.Edges, E{1});
        end
        
        function testFull2D(testCase)
            V = {[0; 1; 2], [0; 1]};
            E = {[1, 2; 2, 3], [1, 2]};
            graph = approx.mesh.SeparableGraph(V, E);
            fullGraph = graph.full();
            
            testCase.verifyClass(fullGraph, 'approx.mesh.Graph');
            testCase.verifyEqual(fullGraph.NVertices, 6);
            testCase.verifyEqual(fullGraph.NEdges, 7);
            testCase.verifyEqual(fullGraph.NDims, 2);
            
            expectedVertices = [
                0, 0;
                1, 0;
                2, 0;
                0, 1;
                1, 1;
                2, 1
            ];
            testCase.verifyEqual(fullGraph.Vertices, expectedVertices);
            
            actualEdges = sort(fullGraph.Edges, 2);
            actualEdges = sortrows(actualEdges);
            expectedEdges = [
                1, 2;
                1, 4;
                2, 3;
                2, 5;
                3, 6;
                4, 5;
                5, 6
            ];
            testCase.verifyEqual(actualEdges, expectedEdges);
        end
        
        function testFull3D(testCase)
            V = {[0; 1], [0; 1], [0; 1]};
            E = {[1, 2], [1, 2], [1, 2]};
            graph = approx.mesh.SeparableGraph(V, E);
            fullGraph = graph.full();
            
            testCase.verifyClass(fullGraph, 'approx.mesh.Graph');
            testCase.verifyEqual(fullGraph.NVertices, 8);
            testCase.verifyEqual(fullGraph.NEdges, 12);
            testCase.verifyEqual(fullGraph.NDims, 3);
        end
        
        function testComplexCase(testCase)
            V = {[0; 1; 2], [0; 1], [0; 1]};
            E = {[1, 2; 2, 3], [1, 2], [1, 2]};
            graph = approx.mesh.SeparableGraph(V, E);
            
            testCase.verifyEqual(graph.NTotalVertices, 12);
            testCase.verifyEqual(graph.NTotalEdges, 20);
            
            fullGraph = graph.full();
            testCase.verifyEqual(fullGraph.NVertices, 12);
            testCase.verifyEqual(fullGraph.NEdges, 20);
        end
    end
end