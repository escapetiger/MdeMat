classdef MacroMicroState < physics.state.State

    properties
        vDisc % Velocity discretization
        nMacroModes % Number of macro modes
        nMicroModes % Number of micro modes
        cache % Precomputed values
        flag % Cache flag
        decomposition % Decomposition pattern
    end

    properties (Dependent)
        nDims % Number of dimensions
        nModes % Number of all modes
        microNodes % Cartesian nodes evaluated by micro modes
        microWeights % Weights of micro modes
        microBasisValues % Basis values of micro modes
        microInflowNodes % Inflow Cartesian nodes evaluate by micro modes
        microOutflowNodes % Outflow Cartesian nodes evaluate by micro modes
        microInflowIndices % Indices of inflow nodes
        microOutflowIndices % Indices of outflow nodes
    end

    methods
        function obj = MacroMicroState(xDisc, vDisc, nMacroModes, nMicroModes)
            obj@physics.state.State(xDisc);
            obj.vDisc = vDisc;
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
            obj.nMacroModes = nMacroModes;
            obj.nMicroModes = nMicroModes;
            obj.clean();
        end

        function n = get.nDims(obj)
            n1 = obj.xDisc.element.geometry.nDims;
            n2 = obj.vDisc.element.geometry.nDims;
            n = min(n1, n2);
        end

        function n = get.nModes(obj)
            n = obj.nMacroModes + obj.nMicroModes;
        end

        function V = get.microNodes(obj)
            d = obj.vDisc.element.geometry.nDims;
            V = obj.vDisc.element.volumeData.nodes;
            if d == 2
                V = obj.vDisc.element.geometry.sphericalToCartesian(V);
            elseif d == 3
                V(1, :) = acos(-V(1, :));
                V = obj.vDisc.element.geometry.sphericalToCartesian(V);
            end
            V = V(1:obj.xDisc.element.nDims, :);
        end

        function w = get.microWeights(obj)
            w = obj.vDisc.element.volumeData.weights;
        end

        function B = get.microBasisValues(obj)
            B = obj.vDisc.element.volumeData.values;
        end

        function V = get.microInflowNodes(obj)
            nd = obj.nDims;
            V = cell(1, 2*nd);
            for d = 1:nd
                VV = obj.microNodes;
                V{2*d-1} = VV(:, VV(d, :) > 0);
                V{2*d} = VV(:, VV(d, :) < 0);
            end
        end

        function V = get.microOutflowNodes(obj)
            nd = obj.nDims;
            V = cell(1, 2*nd);
            for d = 1:nd
                VV = obj.microNodes;
                V{2*d-1} = VV(:, VV(d, :) < 0);
                V{2*d} = VV(:, VV(d, :) > 0);
            end
        end

        function I = get.microInflowIndices(obj)
            nd = obj.nDims;
            I = cell(1, 2*nd);
            for d = 1:nd
                VV = obj.microNodes;
                I{2*d-1} = find(VV(d, :) > 0);
                I{2*d} = find(VV(d, :) < 0);
            end
        end

        function I = get.microOutflowIndices(obj)
            nd = obj.nDims;
            I = cell(1, 2*nd);
            for d = 1:nd
                VV = obj.microNodes;
                I{2*d-1} = find(VV(d, :) < 0);
                I{2*d} = find(VV(d, :) > 0);
            end
        end

        function lazyProject(obj, f, varargin)
            nGlobalDofs = obj.xDisc.space.nGlobalDofs;

            if length(varargin) >= 1 && ~isempty(varargin{1})
                XX = varargin{1};
            else
                I = obj.xDisc.mesh.allElementMultiIndices;
                xRef = obj.xDisc.element.volumeData.nodes;
                XX = obj.xDisc.mesh.collocate(xRef, I);
            end

            if nargin(f) == 1
                obj.cache.C = obj.xDisc.space.project(f, XX);
                obj.cache.C = reshape(obj.cache.C, nGlobalDofs, []);
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
                obj.cache.F = obj.xDisc.space.project(f, X, V);
                obj.cache.F = reshape(obj.cache.F, nGlobalDofs, []);
                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);
                return;
            end

            if length(varargin) >= 3 && ~isempty(varargin{3})
                t = varargin{3};
            else
                t = 0;
            end

            if nargin(f) == 3
                obj.cache.F = obj.xDisc.space.project(f, X, V, t);
                obj.cache.F = reshape(obj.cache.F, nGlobalDofs, []);
                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);
            end
        end

        function lazyFunctionEvaluate(obj, f, varargin)
            if length(varargin) >= 1 && ~isempty(varargin{1})
                XX = varargin{1};
            else
                I = obj.xDisc.mesh.allElementMultiIndices;
                xRef = obj.xDisc.element.volumeData.nodes;
                XX = obj.xDisc.mesh.collocate(xRef, I);
            end

            if nargin(f) == 1
                obj.cache.C = f(XX);
                obj.cache.C = reshape(obj.cache.C, size(XX, 2), []);
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
                obj.cache.F = f(X, V);
                obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);
                return;
            end

            if length(varargin) >= 3 && ~isempty(varargin{3})
                t = varargin{3};
            else
                t = 0;
            end

            if nargin(f) == 3
                obj.cache.F = f(X, V, t);
                obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);
            end
        end

        function lazyTraceEvaluate(obj, i, f, varargin)
            m = obj.nMacroModes;
            n = obj.nMicroModes;

            I = obj.xDisc.mesh.allElementMultiIndices;
            dim = ceil(i / 2);
            if mod(i, 2) == 1
                e = I(:, dim) == 1;
            else
                e = I(:, dim) == obj.xDisc.mesh.resolution(dim);
            end
            I = I(e, :);

            xRef = obj.xDisc.element.fluxData(i).nodes;
            if length(varargin) >= 1 && ~isempty(varargin{1})
                XX = varargin{1};
            else
                XX = obj.xDisc.mesh.collocate(xRef, I);
            end

            if length(varargin) >= 2 && ~isempty(varargin{2})
                VV = varargin{2};
            else
                VV = obj.microNodes;
            end

            X = kron(ones(1, size(VV, 2)), XX);
            V = kron(VV, ones(1, size(XX, 2)));
            if nargin(f) == 2
                nl = obj.xDisc.space.nLocalDofs;
                L = obj.xDisc.space.mesh.indexer.multiToLinear(I);
                K = bsxfun(@plus, (L(:).' - 1)*nl, (1:nl).');
                obj.cache.F = f(X, V);
                obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                if m > 0 && n == 0
                    Jo = obj.microOutflowIndices{i};
                    Vo = obj.vDisc.element.volumeData.nodes(:, Jo);
                    Co = obj.xDisc.space.evaluate(L, xRef, obj.dofs.U(K, :));
                    P = obj.vDisc.element.projector;
                    Uo = P.project(Co.');
                    obj.cache.F(:, Jo) = obj.vDisc.space.evaluate(Vo, Uo, 0);
                elseif m == 0 && n > 0
                    % Nothing to do
                else
                    Jo = obj.microOutflowIndices{i};
                    Vo = obj.vDisc.element.volumeData.nodes(:, Jo);
                    Co = obj.xDisc.space.evaluate(L, xRef, obj.dofs.U(K, :));
                    Go = obj.xDisc.space.evaluate(L, xRef, obj.dofs.G(K, Jo));
                    P = obj.vDisc.element.projector;
                    Uo = P.project(Co.');
                    obj.cache.F(:, Jo) = obj.vDisc.space.evaluate(Vo, Uo, Go);
                end

                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);

                return;
            end

            if length(varargin) >= 3 && ~isempty(varargin{3})
                t = varargin{3};
            else
                t = 0;
            end

            if nargin(f) == 3
                nl = obj.xDisc.space.nLocalDofs;
                L = obj.xDisc.space.mesh.indexer.multiToLinear(I);
                K = bsxfun(@plus, (L(:).' - 1)*nl, (1:nl).');
                obj.cache.F = f(X, V, t);
                obj.cache.F = reshape(obj.cache.F, size(XX, 2), []);
                if m > 0 && n == 0
                    Jo = obj.microOutflowIndices{i};
                    Vo = obj.vDisc.element.volumeData.nodes(:, Jo);
                    Co = obj.xDisc.space.evaluate(L, xRef, obj.dofs.U(K, :));
                    P = obj.vDisc.element.projector;
                    Uo = P.project(Co.');
                    obj.cache.F(:, Jo) = obj.vDisc.space.evaluate(Vo, Uo, 0);
                elseif m == 0 && n > 0
                    % Nothing to do
                else
                    Jo = obj.microOutflowIndices{i};
                    Vo = obj.vDisc.element.volumeData.nodes(:, Jo);
                    Co = obj.xDisc.space.evaluate(L, xRef, obj.dofs.U(K, :));
                    Go = obj.xDisc.space.evaluate(L, xRef, obj.dofs.G(K, Jo));
                    P = obj.vDisc.element.projector;
                    Uo = P.project(Co.');
                    obj.cache.F(:, Jo) = obj.vDisc.space.evaluate(Vo, Uo, Go);
                end
                    
                B = obj.vDisc.element.volumeData.values;
                w = obj.vDisc.element.volumeData.weights;
                S = obj.vDisc.element.projector.embed(obj.cache.F.', B, w);
                obj.cache.C = obj.vDisc.element.projector.project(S);
            end
        end

        function U = macroProject(obj, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b01)
                obj.lazyProject(f, varargin{:});
            end
            M = obj.vDisc.element.projector.mass;
            U = (M * obj.cache.C).';
%             U = U(:, 1:obj.nMacroModes);
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function G = microProject(obj, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b10)
                obj.lazyProject(f, varargin{:});
            end
            v = obj.vDisc.element.volumeData.nodes;
            G = obj.vDisc.space.direction(obj.cache.F, v, obj.cache.C);
%             G = G(:, 1:obj.nMicroModes);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function U = macroFunctionEvaluate(obj, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b01)
                obj.lazyFunctionEvaluate(f, varargin{:});
            end
            M = obj.vDisc.element.projector.mass;
            U = (M * obj.cache.C).';
%             U = U(:, 1:obj.nMacroModes);
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function G = microFunctionEvaluate(obj, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b10)
                obj.lazyFunctionEvaluate(f, varargin{:});
            end
            v = obj.vDisc.element.volumeData.nodes;
            G = obj.vDisc.space.direction(obj.cache.F, v, obj.cache.C);
%             G = G(:, 1:obj.nMicroModes);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function U = macroTraceEvaluate(obj, i, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b01)
                obj.lazyTraceEvaluate(i, f, varargin{:});
            end
            M = obj.vDisc.element.projector.mass;
            U = (M * obj.cache.C).';
            U = U(:, 1:obj.nMacroModes);
            obj.flag = bitor(obj.flag, 0b01);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function G = microTraceEvaluate(obj, i, f, varargin)
            if isempty(obj.cache) || bitand(obj.flag, 0b10)
                obj.lazyTraceEvaluate(i, f, varargin{:});
            end
            v = obj.vDisc.element.volumeData.nodes;
            G = obj.vDisc.space.direction(obj.cache.F, v, obj.cache.C);
            G = G(:, 1:obj.nMicroModes);
            obj.flag = bitor(obj.flag, 0b10);
            if bitand(obj.flag, 0b01) && bitand(obj.flag, 0b10)
                obj.clean(); 
            end
        end

        function C = mean(obj, f)
            C = obj.xDisc.space.average(f);
        end

        function obj = clean(obj)
            obj.cache = [];
            obj.flag = 0b00;
        end
    end
end
