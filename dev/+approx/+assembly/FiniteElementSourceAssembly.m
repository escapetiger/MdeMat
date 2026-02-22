classdef FiniteElementSourceAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTSOURCEASSEMBLY Assembly class for finite element source
    % terms.
    %
    %   FiniteElementSourceAssembly provides functionality for assembling
    %   vector contributions from source terms in finite element
    %   formulations. This class handles the integration of source
    %   functions over finite elements using volume quadrature rules.
    %
    % See also:
    %   approx.assembly.FiniteElementAssembly

    methods
        function G = assembleVector(obj, func, options)
            % ASSEMBLEVECTOR Assemble vector contribution from source
            % terms.
            %
            %   G = assembleVector(obj, func) assembles the finite element
            %   vector from a source function @a func using volume
            %   integration.
            %
            %   G = assembleVector(obj, func, args=args) passes additional
            %   arguments @a args to the function evaluation.

            arguments
                obj approx.assembly.FiniteElementSourceAssembly
                func{mustBeA(func, {'numeric', 'function_handle'})}
                options.args = {}
            end

            nl = obj.Space.NLocalDofs;
            V = obj.Space.Element.Volume.Values;
            W = diag(obj.Space.Element.Volume.Weights);
            nq = obj.Space.Element.Volume.NPoints;
            xRef = obj.Space.Element.Volume.Nodes;
            U = obj.Space.feval(func, xRef, args = options.args);
            U = reshape(U, nq, []);
            L = V * W * U;
            L = obj.Space.Element.Approximator.project(L);
            L = reshape(L, nl, []);
            G = obj.createVector(L);
        end
    end
end
