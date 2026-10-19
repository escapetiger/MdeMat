classdef VolumeAssembly < approx.assembly.Assembly
    % VOLUMEASSEMBLY Assembly for volume terms in finite element methods.
    %
    %   VolumeAssembly handles the assembly of volume integral terms that
    %   arise in finite element discretizations. These terms represent
    %   the weak form contributions from element interiors, such as
    %   diffusion, advection, and reaction terms in partial differential
    %   equations.
    %
    %   The class supports both uniform and non-uniform meshes,
    %   automatically handling the scaling of coefficients based on local
    %   mesh spacing. Volume terms are assembled using element-level
    %   operators and accumulated into global sparse matrices using
    %   appropriate connectivity information.
    %
    %   This assembly process is fundamental to finite element methods, as
    %   volume integrals typically represent the dominant computational
    %   cost and the primary source of system matrix entries in most
    %   discretizations.
    %
    % See also:
    %   approx.assembly.Assembly, approx.assembly.FluxAssembly,
    %   approx.space.MeshSpace, approx.element.ElementOperator

    methods
        function b = scaleConstant(obj, dim, coe)
            % SCALECONSTANT Scale constant coefficients on the mesh.
            %
            %   b = scaleConstant(obj, dim, coe) computes scaled
            %   coefficients on the mesh by dividing input coefficients by
            %   the appropriate mesh spacing in the specified dimension.
            %   This scaling is essential for maintaining proper units and
            %   numerical consistency in finite element formulations.
            %
            %   The scaling accounts for the Jacobian of the transformation
            %   from reference elements to physical elements, ensuring that
            %   the assembled operators have the correct physical
            %   dimensions and magnitude.
            %
            % Inputs:
            %   obj - The VolumeAssembly object
            %   dim - Dimension index for scaling (positive integer)
            %   coe - Input coefficients (scalar or vector)
            %
            % Outputs:
            %   b - Scaled coefficient vector (nElements x 1)

            switch class(obj.space.mesh)
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
            %   dimension dim with the coefficient coe. The method extracts
            %   element-level gradient operators and assembles them into
            %   global sparse matrix format using triplet representation.
            %
            %   The assembled operator corresponds to the weak form of
            %   first-order derivatives, commonly arising in advection
            %   terms or as part of second-order operators like diffusion
            %   when combined with appropriate test function derivatives.
            %
            % Inputs:
            %   obj - The VolumeAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for volume terms

            ne = obj.space.nTotalElements;
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;
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
            % SCALECONSTANTONUNIFORMGRID Scale constant coefficients on
            % uniform grid.

            ne = obj.space.nTotalElements;
            h = obj.space.mesh.spacings{dim};
            b = ones(1, ne) .* coe ./ h;
            b = b(:);
        end

        function b = scaleConstantOnNonuniformGrid(obj, dim, coe)
            % SCALECONSTANTONNONUNIFORMGRID Scale constant coefficients on
            % non-uniform grid.

            h = obj.space.mesh.spacings;
            [h{:}] = ndgrid(h{:});
            b = coe ./ h{dim}(:);
            b = b(:);
        end
    end
end