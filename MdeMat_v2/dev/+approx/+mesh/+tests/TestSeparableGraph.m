classdef TestSeparableGraph < matlab.unittest.TestCase
    
    methods (Test)
        function testConstructor(testCase)
            V1D = {[0; 1; 2; 3]};
            E1D = {[1, 2; 2, 3; 3, 4]};
            graph1D = approx.mesh.SeparableGraph(V1D, E1D);
            testCase.verifyEqual(graph1D.vertices, V1D);
            testCase.verifyEqual(graph1D.edges, E1D);
            
            V2D = {[0; 1; 2], [0; 1; 2]};
            E2D = {[1, 2; 2, 3], [1, 2; 2, 3]};
            graph2D = approx.mesh.SeparableGraph(V2D, E2D);
            testCase.verifyEqual(graph2D.vertices, V2D);
            testCase.verifyEqual(graph2D.edges, E2D);
            
            V3D = {[0; 1], [0; 1], [0; 1]};
            E3D = {[1, 2], [1, 2], [1, 2]};
            graph3D = approx.mesh.SeparableGraph(V3D, E3D);
            testCase.verifyEqual(graph3D.vertices, V3D);
            testCase.verifyEqual(graph3D.edges, E3D);
        end
        
        function testConstructorValidation(testCase)
            V = {[0; 1; 2]};
            E = {[1, 2; 2, 3]};
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(), ...
                'approx:mesh:SeparableGraph:InvalidInput');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V), ...
                'approx:mesh:SeparableGraph:InvalidInput');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph('not cell', E), ...
                'approx:mesh:SeparableGraph:InvalidVertices');
            
            testCase.verifyError(...
                @() approx.mesh.SeparableGraph(V, 'not cell'), ...
                'approx:mesh:SeparableGraph:InvalidEdges');
            
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
            testCase.verifyEqual(graph1D.nDims, 1);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1], [0; 1]}, {[1, 2], [1, 2]});
            testCase.verifyEqual(graph2D.nDims, 2);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.nDims, 3);
        end
        
        function testNVertices(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.nVertices, 4);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.nVertices, [3, 3]);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.nVertices, [2, 2, 2]);
        end
        
        function testNTotalVertices(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.nTotalVertices, 4);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.nTotalVertices, 9);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.nTotalVertices, 8);
        end
        
        function testNEdges(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.nEdges, 3);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.nEdges, [2, 2]);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.nEdges, [1, 1, 1]);
        end
        
        function testNTotalEdges(testCase)
            graph1D = approx.mesh.SeparableGraph({[0; 1; 2; 3]}, {[1, 2; 2, 3; 3, 4]});
            testCase.verifyEqual(graph1D.nTotalEdges, 3);
            
            graph2D = approx.mesh.SeparableGraph({[0; 1; 2], [0; 1; 2]}, {[1, 2; 2, 3], [1, 2; 2, 3]});
            testCase.verifyEqual(graph2D.nTotalEdges, 12);
            
            graph3D = approx.mesh.SeparableGraph({[0; 1], [0; 1], [0; 1]}, {[1, 2], [1, 2], [1, 2]});
            testCase.verifyEqual(graph3D.nTotalEdges, 12);
        end
        
        function testFull1D(testCase)
            V = {[0; 1; 2; 3]};
            E = {[1, 2; 2, 3; 3, 4]};
            graph = approx.mesh.SeparableGraph(V, E);
            fullGraph = graph.full();
            
            testCase.verifyClass(fullGraph, 'approx.mesh.Graph');
            testCase.verifyEqual(fullGraph.nVertices, 4);
            testCase.verifyEqual(fullGraph.nEdges, 3);
            testCase.verifyEqual(fullGraph.vertices, V{1});
            testCase.verifyEqual(fullGraph.edges, E{1});
        end
        
        function testFull2D(testCase)
            V = {[0; 1; 2], [0; 1]};
            E = {[1, 2; 2, 3], [1, 2]};
            graph = approx.mesh.SeparableGraph(V, E);
            fullGraph = graph.full();
            
            testCase.verifyClass(fullGraph, 'approx.mesh.Graph');
            testCase.verifyEqual(fullGraph.nVertices, 6);
            testCase.verifyEqual(fullGraph.nEdges, 7);
            testCase.verifyEqual(fullGraph.nDims, 2);
            
            expectedVertices = [
                0, 0;
                1, 0;
                2, 0;
                0, 1;
                1, 1;
                2, 1
            ];
            testCase.verifyEqual(fullGraph.vertices, expectedVertices);
            
            actualEdges = sort(fullGraph.edges, 2);
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
            testCase.verifyEqual(fullGraph.nVertices, 8);
            testCase.verifyEqual(fullGraph.nEdges, 12);
            testCase.verifyEqual(fullGraph.nDims, 3);
        end
        
        function testComplexCase(testCase)
            V = {[0; 1; 2], [0; 1], [0; 1]};
            E = {[1, 2; 2, 3], [1, 2], [1, 2]};
            graph = approx.mesh.SeparableGraph(V, E);
            
            testCase.verifyEqual(graph.nTotalVertices, 12);
            testCase.verifyEqual(graph.nTotalEdges, 20);
            
            fullGraph = graph.full();
            testCase.verifyEqual(fullGraph.nVertices, 12);
            testCase.verifyEqual(fullGraph.nEdges, 20);
        end
    end
end