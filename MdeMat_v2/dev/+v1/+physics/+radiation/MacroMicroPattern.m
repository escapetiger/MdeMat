classdef MacroMicroPattern < handle
    % MACROMICROPATTERN

    properties
        macroScale
        microScale
        timeScale
        scatteringScale
        absorptionScale
        epsilon
        reduction
        decomposition
        scaling
        streaming
    end

    properties (Dependent)
        scales % Scaling parameters
    end

    methods
        function obj = MacroMicroPattern(varargin)
            % MACROMICROPATTERN

            p = inputParser;
            addParameter(p, 'epsilon', 1);
            addParameter(p, 'macroScale', []);
            addParameter(p, 'microScale', []);
            addParameter(p, 'timeScale', 0);
            addParameter(p, 'scatteringScale', 0);
            addParameter(p, 'absorptionScale', 0);
            addParameter(p, 'reduction', 1);
            parse(p, varargin{:});

            obj.epsilon = p.Results.epsilon;
            obj.macroScale = p.Results.macroScale;
            obj.microScale = p.Results.microScale;
            if isempty(p.Results.timeScale)
                obj.timeScale = 0;
            else
                obj.timeScale = p.Results.timeScale;
            end
            if isempty(p.Results.scatteringScale)
                obj.scatteringScale = 0;
            else
                obj.scatteringScale = p.Results.scatteringScale;
            end
            if isempty(p.Results.absorptionScale)
                obj.absorptionScale = 0;
            else
                obj.absorptionScale = p.Results.absorptionScale;
            end
            obj.reduction = p.Results.reduction;
            obj.decomposition = [];
            obj.streaming = [];
            obj.scaling = [];
        end

        function s = get.scales(obj)
            s = obj.epsilon .^ obj.decomposition;
        end

        function setDecompositionPattern(obj, d, m)
            obj.decomposition = zeros(1, m + 1);
            if d == 1
                obj.decomposition = 0:m;
            elseif d == 2 && obj.reduction == 1
                k = ceil((m-1)/2);
                l = [0, repelem(1:k, 2)];
                obj.decomposition(1:m) = l(1:m);
                obj.decomposition(m+1) = k + 1;
            else
                k = 0;
                j = 0;
                while j < m
                    s = repmat(k, 1, 2*k+1);
                    i = j+1:min(j+2*k+1, m);
                    obj.decomposition(i) = s(1:length(i));
                    j = j + 2 * k + 1;
                    k = k + 1;
                end
                obj.decomposition(m + 1) = k + 1;
            end
        end

        function setFreeStreamingPattern(obj, state)
            TOL = 1e-13;
            m = state.nMacroModes;
            n = state.nMicroModes;
            d = state.nDims;
            W = diag(state.microWeights);
            U = state.microBasisValues.';
            V = state.microNodes;
            B = U ./ diag(U.'*W*U).';
            L = obj.decomposition(1:m);
            J = ones(1, m);
            J(L + 1 < max(L) + 1) = 0;
            J = diag(J);
            P = [eye(m, m), U.' * W; -B, eye(n, n) - B * U.' * W];
            Q = [eye(m, m) - U.' * W * B, -U.' * W; B, eye(n, n)];
            obj.streaming = zeros(m+n, m+n, d);
            for i = 1:d
                C = diag(V(i, :));
                D = [(eye(m,m)-J)*U.'*W*C*B, zeros(m, n); zeros(n, m), C];
                obj.streaming(:, :, i) = P * D * Q;
            end
            obj.streaming(abs(obj.streaming) < TOL) = 0;

%             T = zeros(m+n, m+n, d);
%             I = eye(n, n);
%             for i = 1:d
%                 C = diag(V(i, :));
%                 T(1:m, 1:m, i) = U.' * W * C * B;
%                 T(1:m, m+1:m+n, i) = J * U.' * W * C;
%                 T(m+1:m+n, 1:m, i) = (I - B * U.' * W) * C * B;
%                 T(m+1:m+n, m+1:m+n, i) = (I - B * J * U.' * W) * C;
%             end
        end
    
        function setScalingPattern(obj, state)
            m = state.nMacroModes;
            n = state.nMicroModes;
            if isempty(obj.macroScale)
                scm = obj.decomposition(1:m);
            else
                scm = obj.macroScale;
            end
            if isempty(obj.microScale)
                scn = obj.decomposition(m+1);
            else
                scn = obj.microScale;
            end
            sct = obj.timeScale;
            scs = obj.scatteringScale;
            sca = obj.absorptionScale;
            sc = [scm, repmat(scn, 1, n)];

            S1 = sct + sc(:);
            S2 = repmat(sc, m+n, 1);
            S2(all(obj.streaming == 0, 3)) = inf;
            S3 = inf(m+n, 1);
            if ~isempty(scs), S3(2:end) = scs + sc(2:end); end
            S4 = inf(m+n, 1);
            if ~isempty(sca), S4 = sca + sc(:); end
            S = [S1, S2, S3, S4];
            S = S - min(S, [], 2);
            obj.scaling = (obj.epsilon.^S) .* (S ~= inf);
        end
    end
end
