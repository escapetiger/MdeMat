classdef TestBundledDistMesh < matlab.unittest.TestCase
    methods (TestClassSetup)
        function enableActualExample(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'toolbox', 'examples', 'geometry')));
        end
    end

    methods (Test)
        function circleMeshHasPhysicalGeometryAndPreservesCaller(testCase)
            originalFigures = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() delete(setdiff( ...
                findall(groot, 'Type', 'figure'), originalFigures)));
            sentinel = figure('Visible', 'off');
            axesHandle = axes(sentinel);
            callerLine = plot(axesHandle, [0, 1], [2, 3]);
            originalPath = path;
            originalRng = rng;
            figures = findall(groot, 'Type', 'figure');

            result = ex04_nonuniform_trimesh(ElementSize=0.3, MakePlots=false, Seed=7);
            testCase.verifyEqual(path, originalPath);
            testCase.verifyEqual(rng, originalRng);
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), figures);
            testCase.verifyEqual(get(groot, 'CurrentFigure'), sentinel);
            testCase.verifyEqual(callerLine.YData, [2, 3]);

            points = result.mesh.Vertices;
            triangles = result.mesh.Elements;
            edge1 = points(triangles(:, 2), :)-points(triangles(:, 1), :);
            edge2 = points(triangles(:, 3), :)-points(triangles(:, 1), :);
            areas = (edge1(:, 1).*edge2(:, 2)-edge1(:, 2).*edge2(:, 1))/2;
            testCase.verifyGreaterThan(areas, 0);
            testCase.verifyLessThan(abs(sum(areas)-pi), 0.1);
            testCase.verifyLessThanOrEqual(max(vecnorm(points, 2, 2)), 1+1e-7);
            testCase.verifyGreaterThan(result.quality.minQuality, 0.2);
            testCase.verifyEqual(result.physicalArea, sum(areas), 'AbsTol', 1e-12);
            testCase.verifyEqual(result.refined.NElements, 4*result.mesh.NElements);
            refinedArea = sum(result.refined.computeElementJacobianDeterminants())/2;
            testCase.verifyEqual(refinedArea, sum(areas), 'AbsTol', 1e-12);
            mapped = reshape(result.physicalPoints, 2, 4, []);
            for index = 1:numel(result.targetElements)
                vertices = triangles(result.targetElements(index), :);
                centroid = mean(points(vertices, :), 1).';
                testCase.verifyEqual(mapped(:, 4, index), centroid, 'AbsTol', 1e-13);
            end
        end

        function sourceDistanceHelperIsPureMatlabAndCorrect(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            source = fullfile(root, 'third_party', 'distmesh', 'src');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source));
            testCase.verifyEqual(which('dsegment'), fullfile(source, 'dsegment.m'));
            testCase.verifyEqual(dsegment([1, 0; 1, 1; 3, 0], [0, 0; 2, 0]), ...
                [0; 1; 1], 'AbsTol', 1e-13);
            testCase.verifyEqual(dcircle([0, 0; 1, 0; 2, 0], 0, 0, 1), ...
                [-1; 0; 1], 'AbsTol', 1e-13);
        end
    end
end
