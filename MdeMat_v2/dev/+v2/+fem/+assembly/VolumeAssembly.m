classdef VolumeAssembly < fem.assembly.Assembly
    % VOLUMEASSEMBLY Assembly for volume terms.

    methods
        function b = scaleConstant(obj, dim, coe)
            % SCALECONSTANT Scale constant coefficients on the mesh.
            %
            %   b = scaleConstant(obj, dim, coe) computes scaled
            %   coefficients on the mesh.
            %
            % Inputs:
            %   obj - The VolumeAssembly object
            %   coe - Coefficients
            %
            % Outputs:
            %   b - Scaled coefficient vector

            switch class(obj.context.mesh)
                case 'approx.mesh.UniformGrid'
                    b = obj.scaleConstantOnUniformGrid(dim, coe);
                case 'approx.mesh.NonuniformGrid'
                    b = obj.scaleConstantOnNonuniformGrid(dim, coe);
            end
        end

        function T = assembleVolumePartial(obj, dim, coe)
            % ASSEMBLEVOLUMEPARTIAL Assemble volume terms for partial
            % operator.
            %
            %   T = assembleVolumePartial(obj, dim, coe) creates triplets
            %   for the volume contribution of partial operators along the
            %   dimension @a dim with the coefficient @a coe.
            %
            % Inputs:
            %   obj - The VolumeAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for volume terms

            ne = obj.context.nTotalElements;
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;
            F = G.volumeData(dim).matrix;

            T = zeros(nl^2*ne, 3);
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (0:ne - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(:), v(:));
        end
    end

    methods (Access = private)
        function b = scaleConstantOnUniformGrid(obj, dim, coe)
            % SCALECONSTANTONUNIFORMGRID Scale (piecewise) constant
            % coefficients on the uniform grid.
            ne = obj.context.nTotalElements;
            h = obj.context.mesh.spacings{dim};
            b = ones(1, ne) .* coe ./ h;
            b = b(:);
        end

        function b = scaleConstantOnNonuniformGrid(obj, dim, coe)
            % SCALECONSTANTONNONUNIFORMGRID Scale (piecewise) constant
            % coefficients on the nonuniform grid.
            h = obj.context.mesh.spacings;
            [h{:}] = ndgrid(h{:});
            b = coe ./ h{dim}(:);
            b = b(:);
        end
    end
end
