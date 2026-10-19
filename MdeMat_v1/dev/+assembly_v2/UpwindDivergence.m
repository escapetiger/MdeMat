classdef UpwindDivergence < fem.assembly.Operator
    % DIVEREGNCE

    properties
        partial
    end

    methods
        function obj = UpwindDivergence(volume, flux)
            obj.partial = fem.assembly.L2Partial(volume, flux);
        end

        function A = linear(obj, c, a)
            if nargin < 3 || isempty(a), a = c; end

            n = obj.nGlobalDofs;
            A = sparse(n, n);
            for d = 1:obj.nDims
                if c(d) == 0, continue; end
                if a(d) > 0
                    T = obj.partial.assemble(d, c(d), 1);
                else
                    T = obj.partial.assemble(d, c(d), 2);
                end
                A = A + T;
            end
        end
    
        function b = linearBc(obj, c, a, f)
            n = obj.partial.volume.nGlobalDofs;
            b = sparse(n, 1);
            for d = 1:obj.partial.volume.nDims
                if c(d) == 0, continue; end
                b = b + obj.partial.assembleBc(d, c(d), a(d), f);
            end
        end

        function b = linearPartialBc(obj, d, c, a)

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{d};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{d}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                if a > 0
                    T = obj.assembleLeftTrace(d, c);
                else
                    T = obj.assembleRightTrace(d, c);
                end
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            b = obj.sparseFromTriplets(T, n, 1);
        end
    end
end
