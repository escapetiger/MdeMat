classdef GridVolumeAssembly < fem.assembly.Assembly
    %VOLUMEASSEMBLY
    
    methods
        function b = assembleConstant(obj, c)
            switch class(obj.mesh)
                case 'approx.mesh.UniformGrid'
                    b = obj.assembleConstantOnUniformGrid(c);
                case 'approx.mesh.NonuniformGrid'
                    b = obj.assembleConstantOnNonuniformGrid(c);
            end
        end

        function b = assembleConstantOnBoundary(obj, c)
        end
        
        function T = assembleVolumePartial(obj, d, c)
            ne = obj.nTotalElements;
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;
            F = G.volumeData(d).matrix;

            T = zeros(np^2*ne, 3);
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (0:ne - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (0:ne - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(:), v(:));
        end
    end

    methods
        function b = assembleConstantOnUniformGrid(obj, c)
            ne = obj.nTotalElements;
            h = obj.mesh.spacings{d};
            b = ones(1, ne) .* c ./ h;
            b = b(:);
        end

        function b = assembleConstantOnNonuniformGrid(obj, c)
            h = obj.mesh.spacings;
            [h{:}] = ndgrid(h{:});
            b = c ./ h{d}(:);
            b = b(:);
        end
    end
end

