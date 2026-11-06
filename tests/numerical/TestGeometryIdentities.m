classdef TestGeometryIdentities < matlab.unittest.TestCase
    methods (Test)
        function planeOffsetAndNormalScaling(testCase)
            intersector = mdemat.common.geometry.SegmentByPlaneIntersector(2);
            start = [0, 0; 0, 1];
            finish = [3, 4; 0, 1];
            normals = [1, 3; 0, 0];
            offsets = [-1, -6];
            [points, found] = intersector.intersect(start, finish, normals, offsets);
            testCase.verifyTrue(all(found));
            testCase.verifyEqual(points, [1, 2; 0, 1], 'AbsTol', 2e-14);
            [scaled, found] = intersector.intersect(start, finish, ...
                1e-5*normals, 1e-5*offsets);
            testCase.verifyTrue(all(found));
            testCase.verifyEqual(scaled, points, 'AbsTol', 2e-14);
        end

        function sharedEndpointDoesNotHideCollinearOverlap(testCase)
            intersector = mdemat.common.geometry.SegmentBySegmentIntersector(2);
            for pair = {[0, 1; 0, 0], [1, 2; 0, 0]}
                testCase.verifyError(@() intersector.intersect( ...
                    [0; 0], [2; 0], pair{1}(:, 1), pair{1}(:, 2)), ...
                    'mdemat:common:geometry:SegmentBySegmentIntersector:CollinearOverlap');
            end
            [points, found] = intersector.intersect([0; 0], [1; 0], [2; 0], [1; 0]);
            testCase.verifyTrue(found);
            testCase.verifyEqual(points, [1; 0]);
        end

        function nonuniformMapsAndRefinementPreservePhysicalMeasure(testCase)
            mesh = mdemat.data.mesh.Grid({[-2; -0.3; 1], [0, 0.2, 2, 3]});
            jacobians = mesh.computeElementJacobians();
            inverses = mesh.computeElementInverseJacobians();
            determinants = mesh.computeElementJacobianDeterminants();
            testCase.verifyEqual(sum(determinants), 9, 'AbsTol', 3e-14);
            for element = 1:mesh.NElements
                mapped = mesh.collocate([-0.5, 0.5; -0.5, 0.5], element);
                testCase.verifyEqual(mapped(:, 2)-mapped(:, 1), ...
                    diag(jacobians(:, :, element)), 'AbsTol', 2e-14);
                testCase.verifyEqual(jacobians(:, :, element)*inverses(:, :, element), ...
                    eye(2), 'AbsTol', 2e-14);
            end
            fine = mesh.refine(2);
            testCase.verifyEqual(sum(fine.computeElementJacobianDeterminants()), ...
                sum(determinants), 'AbsTol', 3e-14);
            for invalid = {{[0, 1, 1]}, {[1, 0]}, {[0, NaN]}, {0}}
                testCase.verifyError(@() mdemat.data.mesh.Grid(invalid{1}), ...
                    'MdeMat:Grid:InvalidVertices');
            end
        end

        function sphericalCoordinatesPreserveTranslationAndRadius(testCase)
            for dimension = 2:5
                center = (1:dimension)/7;
                sphere = mdemat.data.geometry.Sphere(center, 2.3);
                angles = repmat(linspace(0.2, 2.8, 9), dimension-1, 1);
                angles(end, :) = linspace(0.1, 6.1, 9);
                points = sphere.sphericalToCartesian(angles);
                testCase.verifyEqual(vecnorm(points-center(:), 2, 1), ...
                    2.3*ones(1, 9), 'AbsTol', 2e-14);
                reconstructed = sphere.sphericalToCartesian(sphere.cartesianToSpherical(points));
                testCase.verifyEqual(reconstructed, points, 'AbsTol', 2e-13);
            end
        end
    end
end
