classdef MacroMicroState < physics.State

    properties
        cache % Precomputed values
        flag
    end

    properties (Dependent)
        nDims % Number of dimensions
        nMacroModes % Number of macro modes
        nMicroModes % Number of micro modes
        nModes % Number of all modes
        microNodes % Cartesian nodes evaluated by micro modes 
        microWeights % Weights of micro modes
        microBasisValues % Basis values of micro modes
    end

    methods
        function obj = MacroMicroState(disc)
            obj@physics.State(disc);
            obj.dofs = struct( ...
                'U', [], ... %< Macro coefficients
                'G', [] ... %< Micro coefficients
                );
            obj.coefs = struct( ...
                'CS', [], ... %< Scattering coefficients
                'CA', [], ... %< Absorption coefficients
                'QU', [], ... %< Macro source coefficients
                'QG', [] ... %< Micro source coefficients
                );
            obj.clean();
        end

        function n = get.nDims(obj)
            n1 = obj.disc.x.fe.geometry.nDims;
            n2 = obj.disc.v.fe.geometry.nDims;
            n = min(n1, n2);
        end

        function n = get.nMacroModes(obj)
            n = size(obj.dofs.U, 2);
        end

        function n = get.nMicroModes(obj)
            n = size(obj.dofs.G, 2);
        end

        function n = get.nModes(obj)
            n = obj.nMacroModes + obj.nMicroModes;
        end

        function V = get.microNodes(obj)
            d = obj.disc.v.fe.geometry.nDims;
            V = obj.disc.v.fe.volumeData.nodes;
            if d == 2
                V = obj.disc.v.fe.geometry.sphericalToCartesian(V);
            elseif d == 3
                V(1, :) = acos(-V(1, :));
                V = obj.disc.v.fe.geometry.sphericalToCartesian(V);
            end
            V = V(1:obj.disc.x.fe.nDims, :);
        end

        function w = get.microWeights(obj)
            w = obj.disc.v.fe.volumeData.weights;
        end

        function B = get.microBasisValues(obj)
            B = obj.disc.v.fe.volumeData.values;
        end

        function precompute(obj, mode, f, varargin)
            nGlobalDofs = obj.disc.x.space.nDofs;
            if length(varargin) >= 1 && ~isempty(varargin{1})
                XX = varargin{1};
            else
                xRef = obj.disc.x.fe.volumeData.nodes;
                I = obj.disc.x.mesh.allElementMultiIndices;
                XX = obj.disc.x.mesh.collocate(xRef, I);
            end

            if nargin(f) == 1
                if mode == 0
                    obj.cache.C = obj.disc.x.space.project(f, XX);
                    obj.cache.C = reshape(obj.cache.C, nGlobalDofs, []);
                else
                    obj.cache.C = f(XX);
                    obj.cache.C = reshape(obj.cache.C, size(XX, 2), []);
                end
                return;
            end

            if length(varargin) >= 2 && ~isempty(varargin{2})
                VV = varargin{2};
            else
                VV = obj.microNodes;
            end

            X = kron(ones(1, size(VV, 2)), XX);
            V = kron(VV, ones(1, size(XX, 2)));
            if nargin(f) == 2
                if mode == 0
                    obj.cache.F = obj.disc.x.space.project(f, X, V);
                    obj.cache.F = reshape(obj.cache.F, nGlobalDofs, []);
                else
                    obj.cache.F = f(X, V);
                    obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                end
                B = obj.disc.v.fe.volumeData.values;
                w = obj.disc.v.fe.volumeData.weights;
                S = obj.disc.v.fe.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.disc.v.fe.projector.project(S);
                return;
            end
            
            if length(varargin) >= 3 && ~isempty(varargin{3})
                t = varargin{3};
            else
                t = 0;
            end
            if nargin(f) == 3
                if mode == 0
                    obj.cache.F = obj.disc.x.space.project(f, X, V, t);
                    obj.cache.F = reshape(obj.cache.F, nGlobalDofs, []);
                else
                    obj.cache.F = f(X, V, t);
                    obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                end
                B = obj.disc.v.fe.volumeData.values;
                w = obj.disc.v.fe.volumeData.weights;
                S = obj.disc.v.fe.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.disc.v.fe.projector.project(S);
            end
        end
        
        function precomputeTrace(obj, mode, i, f, varargin)
            
            %< Inflow

            %< Outflow
        end

        function F = kineticProject(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(0, f, varargin{:});
            end
            F = obj.cache.F;
        end

        function U = macroProject(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(0, f, varargin{:});
            end
            M = obj.disc.v.fe.projector.mass;
            U = (M * obj.cache.C).';
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function G = microProject(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(0, f, varargin{:});
            end
            v = obj.disc.v.fe.volumeData.nodes;
            G = obj.disc.v.space.direction(obj.cache.F, v, obj.cache.C);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function F = kineticEvaluate(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(1, f, varargin{:});
            end
            F = obj.cache.F;
        end

        function U = macroEvaluate(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(1, f, varargin{:});
            end
            M = obj.disc.v.fe.projector.mass;
            U = (M * obj.cache.C).';
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function G = microEvaluate(obj, f, varargin)
            if isempty(obj.cache)
                obj.precompute(1, f, varargin{:});
            end
            v = obj.disc.v.fe.volumeData.nodes;
            G = obj.disc.v.space.direction(obj.cache.F, v, obj.cache.C);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function U = macroTraceProject(obj, i, f, varargin)
            d = ceil(i / 2);
            if isempty(obj.cache)
                obj.precomputeTrace(0, i, f, varargin{:});
            end
            M = obj.disc.v.fe.projector.mass;
            U = (M * obj.cache.C).';
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function G = microTraceProject(obj, i, f, varargin)
            if isempty(obj.cache)
                obj.precompute(0, f, varargin{:});
            end
            v = obj.disc.v.fe.volumeData.nodes;
            G = obj.disc.v.space.direction(obj.cache.F, v, obj.cache.C);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b11), obj.clean(); end
        end

        function C = mean(obj, f)
            C = obj.disc.x.space.average(f);
        end

        function clean(obj)
            obj.cache = [];
            obj.flag = 0b00;
        end
    
        function [V, I] = getInflowPoints(obj, i)
            d = ceil(i/2);
            V = obj.disc.v.fe.volumdData.nodes;
            if mod(i, 2) == 1
                I = V(d, :) > 0;
                V = V(:, I);
                I = find(I);
            else
                I = V(d, :) < 0;
                V = V(:, I);
                I = find(I);
            end
        end

        function [V, I] = getOutflowPoints(obj, i)
            d = ceil(i/2);
            V = obj.disc.v.fe.volumdData.nodes;
            if mod(i, 2) == 1
                I = V(d, :) < 0;
                V = V(:, I);
                I = find(I);
            else
                I = V(d, :) > 0;
                V = V(:, I);
                I = find(I);
            end
        end
    end
end
