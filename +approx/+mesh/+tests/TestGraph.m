classdef TestGraph < matlab.unittest.TestCase
    
    methods (Test)
        function testEmptyConstructor(testCase)
            graph = approx.mesh.Graph();
            testCase.verifyEmpty(graph.vertices);
            testCase.verifyEmpty(graph.edges);
            testCase.verifyEqual(graph.nDims, 0);
            testCase.verifyEqual(graph.nVertices, 0);
            testCase.verifyEqual(graph.nEdges, 0);
        end
        
        function testConstructor(testCase)
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            edges = [1, 2; 1, 3; 2, 4; 3, 4];
            
            graph = approx.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.vertices, vertices);
            testCase.verifyEqual(graph.edges, edges);
        end
        
        function testSetVertices(testCase)
            graph = approx.mesh.Graph();
            
            vertices1D = [0; 1; 2];
            graph.setVertices(vertices1D);
            testCase.verifyEqual(graph.vertices, vertices1D);
            
            vertices2D = [0, 0; 1, 0; 0, 1];
            graph.setVertices(vertices2D);
            testCase.verifyEqual(graph.vertices, vertices2D);
            
            vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0];
            graph.setVertices(vertices3D);
            testCase.verifyEqual(graph.vertices, vertices3D);
        end
        
        function testSetVerticesValidation(testCase)
            graph = approx.mesh.Graph();
            
            testCase.verifyError(...
                @() graph.setVertices('not a matrix'), ...
                'approx:mesh:Graph:InvalidVertices');
        end
        
        function testSetEdges(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            graph = approx.mesh.Graph(vertices, []);
            
            edges = [1, 2; 2, 3];
            graph.setEdges(edges);
            testCase.verifyEqual(graph.edges, edges);
        end
        
        function testSetEdgesValidation(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            graph = approx.mesh.Graph(vertices, []);
            
            testCase.verifyError(...
                @() graph.setEdges('not a matrix'), ...
                'approx:mesh:Graph:InvalidEdges');
            
            testCase.verifyError(...
                @() graph.setEdges([1, 2, 3]), ...
                'approx:mesh:Graph:InvalidEdges');
            
            testCase.verifyError(...
                @() graph.setEdges([1, 4]), ...
                'approx:mesh:Graph:InvalidEdgeIndices');
            
            testCase.verifyError(...
                @() graph.setEdges([0, 2]), ...
                'approx:mesh:Graph:InvalidEdgeIndices');
        end
        
        function testNDims(testCase)
            vertices1D = [0; 1; 2];
            edges1D = [1, 2; 2, 3];
            graph1D = approx.mesh.Graph(vertices1D, edges1D);
            testCase.verifyEqual(graph1D.nDims, 1);
            
            vertices2D = [0, 0; 1, 0; 0, 1];
            edges2D = [1, 2; 2, 3];
            graph2D = approx.mesh.Graph(vertices2D, edges2D);
            testCase.verifyEqual(graph2D.nDims, 2);
            
            vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0];
            edges3D = [1, 2; 2, 3];
            graph3D = approx.mesh.Graph(vertices3D, edges3D);
            testCase.verifyEqual(graph3D.nDims, 3);
        end
        
        function testNVertices(testCase)
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            edges = [1, 2; 1, 3];
            graph = approx.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.nVertices, 4);
            
            emptyGraph = approx.mesh.Graph();
            testCase.verifyEqual(emptyGraph.nVertices, 0);
        end
        
        function testNEdges(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            edges = [1, 2; 2, 3];
            graph = approx.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.nEdges, 2);
            
            emptyGraph = approx.mesh.Graph();
            testCase.verifyEqual(emptyGraph.nEdges, 0);
        end
        
        function testComplexGraph(testCase)
            vertices = [0, 0; 1, 0; 1, 1; 0, 1];
            edges = [1, 2; 2, 3; 3, 4; 4, 1];
            
            graph = approx.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.nVertices, 4);
            testCase.verifyEqual(graph.nEdges, 4);
            testCase.verifyEqual(graph.nDims, 2);
            testCase.verifyEqual(graph.vertices, vertices);
            testCase.verifyEqual(graph.edges, edges);
        end
    end
end