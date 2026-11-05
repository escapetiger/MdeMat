classdef FiniteElementSemiLagrangianAssembly < approx.assembly.FiniteElementAssembly
    % SEMILAGRANGIANASSEMBLY

    properties (Constant)
        CLIP_FOR_MATIX = 0b10 % Clipped when assembling matrix
        CLIP_FOR_VECTOR = 0b01 % Clipped when assembling vector
    end

    properties
        bcType % Boundary condition type
        status % Status mask
    end

    methods
        function obj = FiniteElementSemiLagrangianAssembly(space, bcType)
            obj@approx.assembly.FiniteElementAssembly(space);
            obj.setBcType(bcType);
            obj.reset();
        end

        function obj = setBcType(obj, bcType)
            obj.bcType = bcType;
            obj.checkValidBcType();
        end

        function obj = reset(obj)
            obj.status = 0;
        end

        function G = assembleMatrix(obj, field, tBegin, tEnd)
            d = obj.space.nDims;

            %< Parse advection field
            field = reshape(field, d, []);

            %< Clipping
            ht = tBegin - tEnd;
            hx = cellfun(@(h) h(1), obj.space.mesh.spacings);
            obj.space.element.clipVolume(field, ht, hx);

            %< Update status
            obj.status = bitor(obj.status, obj.CLIP_FOR_MATIX);

            %< Assemble global mass matrix
            switch obj.bcType
                case 'periodic'
                    G = obj.assembleMatrixPeriodic();
                case 'dirichlet'
                    G = obj.assembleMatrixDirichlet();
            end
        end


        function G = assembleVector(obj, field, tBegin, tEnd, func)
            if bitand(obj.status, obj.CLIP_FOR_MATIX)
                d = obj.space.nDims;

                %< Parse advection field
                field = reshape(field, d, []);

                %< Clipping
                ht = tBegin - tEnd;
                hx = cellfun(@(h) h(1), obj.space.mesh.spacings);
                obj.space.element.clipVolume(field, ht, hx);

                %< Update status
                obj.status = bitor(obj.status, obj.CLIP_FOR_MATIX);
            end

            %< Skip for periodic boundary condition
            if strcmpi(obj.bcType, 'periodic'), return; end

            G = obj.assembleVectorDirichlet(field, tBegin, tEnd, func);
        end
    end

    methods (Access = protected)
        function tf = checkValidSpace(obj)
            checkValidSpace@approx.assembly.Assembly(obj);
            tf1 = isa(obj.space.element, 'approx.element.L2SLOrthotopeElement');
            tf2 = isa(obj.space.mesh, 'approx.mesh.UniformGrid');
            tf = tf1 & tf2;
            core.except.assert(tf, 'InvalidInput', ...
                'SL assembly requires L2 SL orthotope element and uniform grid.');
        end

        function tf = checkValidBcType(obj)
            tf = ismember(obj.bcType, {'periodic', 'dirichlet'});
            core.except.assert(tf, 'InvalidInput', 'BC type is not supported.');
        end
    end

    methods (Access = private)
        function G = assembleMatrixPeriodic(obj)
            G = obj.createZeros();
            for i = 1:obj.space.element.nVolumePieces
                %< Compute local mass matrix
                W = diag(obj.space.element.volumePieces.eulerian(i).weights); % nq x nq
                V = obj.space.element.volumePieces.eulerian(i).values; % nl x nq
                U = obj.space.element.volumePieces.upstream(i).values; % nl x nq
                L = V * W * U.'; % nl x nl

                %< Add contribution to global mass matrix
                IE = 1:obj.space.nMeshElements;
                ME = obj.space.mesh.indexer.linearToMulti(IE);
                ME = ME + obj.space.element.volumePieces.indexShifts(:, i).';
                JE = obj.space.mesh.indexer.multiToLinear(ME, 'wrap');
                GP = obj.createMatrix(IE, [], JE, [], L);
                G = G + GP;
            end
        end

        function G = assembleMatrixDirichlet(obj)
            d = obj.space.nDims;
            G = obj.createZeros();
            for i = 1:obj.space.element.nVolumePieces
                %< Compute local mass matrix
                W = diag(obj.space.element.volumePieces.eulerian(i).weights); % nq x nq
                V = obj.space.element.volumePieces.eulerian(i).values; % nl x nq
                U = obj.space.element.volumePieces.upstream(i).values; % nl x nq
                L = V * W * U.'; % nl x nl

                %< Find elements that have at least one interior upstream root
                ME = obj.space.mesh.indexer.generate();
                SE = obj.space.element.volumePieces.indexShifts(:, i).';
                MU = ME + SE;
                mask1 = (MU < 1) .* (2 * (1:d) - 1);
                mask2 = (MU > obj.space.mesh.resolution) .* (2 * (1:d));
                mask = mask1 + mask2;
                mask = sum(mask ~= 0, 2) == 0;
                ME = ME(mask, :);

                %< Add contribution to global mass matrix
                IE = obj.space.mesh.indexer.multiToLinear(ME, 'empty');
                JE = obj.space.mesh.indexer.multiToLinear(ME+SE, 'empty');
                GP = obj.createMatrix(IE, [], JE, [], L);
                G = G + GP;
            end
        end

        function G = assembleVectorDirichlet(obj, field, tBegin, tEnd, func)
            np = obj.space.element.nVolumePieces;
            nd = obj.space.nDims;
            ht = tBegin - tEnd;

            G = obj.createZeros(1);
            for i = 1:np
                %< Find elements whose pth volume piece exits domain
                ME = obj.space.mesh.indexer.generate();
                SE = obj.space.element.volumePieces.indexShifts(:, i).';
                MU = ME + SE;
                mask1 = (MU < 1) .* (2 * (1:nd) - 1);
                mask2 = (MU > obj.space.mesh.resolution) .* (2 * (1:nd));
                mask = mask1 + mask2;
                mask = sum(mask ~= 0, 2) > 0;
                ME = ME(mask, :);
                IE = obj.space.mesh.indexer.multiToLinear(ME, 'empty');
                ne = size(ME, 1);

                if ne == 0, return; end

                %< Set up characteristic end points
                nq = obj.space.element.volumePieces.eulerian(i).nPoints;
                X1 = zeros(nd+1, ne*nq);
                X2 = zeros(nd+1, ne*nq);
                Xq = obj.space.element.volumePieces.eulerian(i).nodes;
                X2(1:nd, :) = obj.space.mesh.collocate(Xq, IE);
                X2(nd+1, :) = tBegin;
                X1(1:nd, :) = X2(1:nd, :) - ht * field(:);
                X1(nd+1, :) = tEnd;

                %< Find boundary intersections of characteristics
                intersector = core.geometry.LineSegmentByPlaneIntersector(nd+1);
                bbox = reshape(obj.space.mesh.bbox, 2, []);
                bbox = -bbox(2*(1:nd)-(field > 0));
                Xb = zeros(nd+1, nq*ne);
                for d = 1:nd
                    a = zeros(nd+1, 1);
                    a(d) = 1;
                    [X, TF] = intersector.intersect(X1, X2, a, bbox(d));
                    Xb(:, TF) = X;
                end

                %< Evaluate boundary function and project to element space
                V = obj.space.element.volumePieces.eulerian(i).values; % nl x nq
                W = diag(obj.space.element.volumePieces.eulerian(i).weights); % nq x nq
                U = func(Xb(1:nd, :), Xb(nd+1, :)); % 1 x nq*ne
                U = reshape(U, nq, ne);
                L = V * W * U; % nl x ne
                GP = obj.createVector(IE, [], L);
                G = G + GP;
            end
        end
    end
end
