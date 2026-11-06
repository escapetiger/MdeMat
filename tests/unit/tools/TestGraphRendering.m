classdef TestGraphRendering < matlab.unittest.TestCase
    methods (Test)
        function graphsKeepCoordinatesAndCallerState(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f));
            ax = axes(f);
            for nd = 1:3
                cla(ax); vertices = [zeros(1, nd); ones(1, nd); 2 * ones(1, nd)];
                graph = mdemat.data.mesh.Graph(vertices, [1, 2; 2, 3]);
                hold(ax, 'off'); h = plot(graph, axes=ax, edgeColor='r');
                test.verifyFalse(ishold(ax));
                test.verifyEqual(graph.Vertices, vertices);
                actual = h.Edges.XData; actual = actual(isfinite(actual));
                test.verifyEqual(actual(:), [0; 1; 1; 2]);
                test.verifyEqual(h.Vertices.XData(:), [0; 1; 2]);
                hold(ax, 'on'); plot(graph, axes=ax, showEdges=false);
                test.verifyTrue(ishold(ax));
            end
            bad = mdemat.data.mesh.Graph([zeros(1, 4); ones(1, 4)], [1, 2]);
            test.verifyError(@() plot(bad, axes=ax), 'mdemat:data:mesh:Graph:UnsupportedDimension');
        end
        function separableSamplingPreservesTopologyAndBudget(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            graph = separable(2, 20);
            h = plot(graph, axes=ax, maxVertices=12, showEdges=false);
            test.verifyLessThanOrEqual(numel(h.Vertices.XData), 12);
            test.verifyEqual([min(h.Vertices.XData), max(h.Vertices.XData)], [0, 1]);
            test.verifyEqual(graph.NTotalVertices, 400);
            isolated = mdemat.data.mesh.SeparableGraph({1:200}, {zeros(0, 2)});
            h = plot(isolated, axes=ax, maxVertices=12, showEdges=false);
            test.verifyEqual(numel(h.Vertices.XData), 12);
            test.verifyEmpty(h.Edges);
            test.verifyError(@() plot(separable(2, 2), axes=ax, maxVertices=3), ...
                'mdemat:data:mesh:SeparableGraph:VertexBudgetTooSmall');
            test.verifyError(@() plot(graph, axes=ax, projectionDims=[1, 1]), ...
                'mdemat:data:mesh:SeparableGraph:InvalidProjection');
            h = plot(separable(4, 3), axes=ax, projectionDims=[2, 4], showEdges=false);
            test.verifyEqual(numel(h.Vertices.XData), 9);
        end
        function meshMethodIsInheritedAndPreservesGeometry(test)
            f = figure('Visible', 'off'); other = figure('Visible', 'off');
            guard = onCleanup(@() delete([f, other])); ax = axes(f);
            for nd = 1:3
                mesh = mdemat.data.mesh.Triangulation([zeros(1, nd); eye(nd)], 1:nd+1);
                original = mesh.Vertices;
                for fast = [false, true]
                    cla(ax); hold(ax, 'off');
                    h = plot(mesh, axes=ax, useFastRendering=fast, ...
                        showFaces=true, faceColor='c', showElements=true, elementColor='b');
                    test.verifyFalse(ishold(ax));
                    test.verifyTrue(isgraphics(other));
                    test.verifyEqual(mesh.Vertices, original);
                    test.verifyNotEmpty(h.Vertices);
                    test.verifyNotEmpty(h.Elements);
                    test.verifyTrue(contains(ax.Title.String, sprintf('%dD', nd)));
                    if nd >= 2, test.verifyEqual(h.Elements.Vertices, original); end
                end
            end
            mesh.Elements = [1, 2, 3, 4, 1];
            test.verifyError(@() plot(mesh, axes=ax), ...
                'mdemat:data:mesh:PolytopalMesh:UnsupportedTopology');
        end
        function highDimensionalMeshProjectionPreservesPrincipalScales(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            vertices = [zeros(1, 4); eye(4)];
            mesh = mdemat.data.mesh.Triangulation(vertices, 1:5);
            test.verifyWarning(@() plot(mesh, axes=ax), ...
                'mdemat:data:mesh:PolytopalMesh:DimensionReduction');
            points = findall(ax, 'Type', 'scatter');
            projected = [points.XData(:), points.YData(:), points.ZData(:)];
            scales = svd(vertices - mean(vertices, 1));
            test.verifyEqual(svd(projected), scales(1:3), AbsTol=3e-15);
            test.verifyEqual(mesh.Vertices, vertices);
        end
        function explicitAxesPreserveUnrelatedContent(test)
            f = figure('Visible', 'off'); other = figure('Visible', 'off');
            guard = onCleanup(@() delete([f, other])); ax = axes(f);
            retained = line(ax, [0, 1], [3, 4]); set(groot, 'CurrentFigure', other);
            plot(mdemat.data.mesh.Graph([0; 1], [1, 2]), axes=ax);
            test.verifyTrue(isgraphics(retained));
            test.verifyEqual(get(groot, 'CurrentFigure'), other);
        end
    end
end
function graph = separable(nd, n)
vertices = repmat({linspace(0, 1, n)}, 1, nd);
edges = repmat({[(1:n-1)', (2:n)']}, 1, nd);
graph = mdemat.data.mesh.SeparableGraph(vertices, edges);
end
