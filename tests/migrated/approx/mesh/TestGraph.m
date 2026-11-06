classdef TestGraph < matlab.unittest.TestCase
    
    methods (Test)
        function testEmptyConstructor(testCase)
            graph = mdemat.data.mesh.Graph();
            testCase.verifyEmpty(graph.Vertices);
            testCase.verifyEmpty(graph.Edges);
            testCase.verifyEqual(graph.NDims, 0);
            testCase.verifyEqual(graph.NVertices, 0);
            testCase.verifyEqual(graph.NEdges, 0);
        end
        
        function testConstructor(testCase)
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            edges = [1, 2; 1, 3; 2, 4; 3, 4];
            
            graph = mdemat.data.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.Vertices, vertices);
            testCase.verifyEqual(graph.Edges, edges);
        end
        
        function testSetVertices(testCase)
            graph = mdemat.data.mesh.Graph();
            
            vertices1D = [0; 1; 2];
            graph.setVertices(vertices1D);
            testCase.verifyEqual(graph.Vertices, vertices1D);
            
            vertices2D = [0, 0; 1, 0; 0, 1];
            graph.setVertices(vertices2D);
            testCase.verifyEqual(graph.Vertices, vertices2D);
            
            vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0];
            graph.setVertices(vertices3D);
            testCase.verifyEqual(graph.Vertices, vertices3D);
        end
        
        function testSetVerticesValidation(testCase)
            graph = mdemat.data.mesh.Graph();
            
            testCase.verifyError(...
                @() graph.setVertices('not a matrix'), ...
                'MATLAB:validators:mustBeNumeric');
        end
        
        function testSetEdges(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            graph = mdemat.data.mesh.Graph(vertices, []);
            
            edges = [1, 2; 2, 3];
            graph.setEdges(edges);
            testCase.verifyEqual(graph.Edges, edges);
        end
        
        function testSetEdgesValidation(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            graph = mdemat.data.mesh.Graph(vertices, []);
            
            testCase.verifyError(...
                @() graph.setEdges('not a matrix'), ...
                'MATLAB:validators:mustBeNumeric');
            
            testCase.verifyError(...
                @() graph.setEdges([1, 2, 3]), ...
                'mdemat:data:mesh:Graph:InvalidEdges');
            
            testCase.verifyError(...
                @() graph.setEdges([1, 4]), ...
                'mdemat:data:mesh:Graph:InvalidEdgeIndices');
            
            testCase.verifyError(...
                @() graph.setEdges([0, 2]), ...
                'mdemat:data:mesh:Graph:InvalidEdgeIndices');
        end
        
        function testNDims(testCase)
            vertices1D = [0; 1; 2];
            edges1D = [1, 2; 2, 3];
            graph1D = mdemat.data.mesh.Graph(vertices1D, edges1D);
            testCase.verifyEqual(graph1D.NDims, 1);
            
            vertices2D = [0, 0; 1, 0; 0, 1];
            edges2D = [1, 2; 2, 3];
            graph2D = mdemat.data.mesh.Graph(vertices2D, edges2D);
            testCase.verifyEqual(graph2D.NDims, 2);
            
            vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0];
            edges3D = [1, 2; 2, 3];
            graph3D = mdemat.data.mesh.Graph(vertices3D, edges3D);
            testCase.verifyEqual(graph3D.NDims, 3);
        end
        
        function testNVertices(testCase)
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            edges = [1, 2; 1, 3];
            graph = mdemat.data.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.NVertices, 4);
            
            emptyGraph = mdemat.data.mesh.Graph();
            testCase.verifyEqual(emptyGraph.NVertices, 0);
        end
        
        function testNEdges(testCase)
            vertices = [0, 0; 1, 0; 0, 1];
            edges = [1, 2; 2, 3];
            graph = mdemat.data.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.NEdges, 2);
            
            emptyGraph = mdemat.data.mesh.Graph();
            testCase.verifyEqual(emptyGraph.NEdges, 0);
        end
        
        function testComplexGraph(testCase)
            vertices = [0, 0; 1, 0; 1, 1; 0, 1];
            edges = [1, 2; 2, 3; 3, 4; 4, 1];
            
            graph = mdemat.data.mesh.Graph(vertices, edges);
            testCase.verifyEqual(graph.NVertices, 4);
            testCase.verifyEqual(graph.NEdges, 4);
            testCase.verifyEqual(graph.NDims, 2);
            testCase.verifyEqual(graph.Vertices, vertices);
            testCase.verifyEqual(graph.Edges, edges);
        end
    end
end