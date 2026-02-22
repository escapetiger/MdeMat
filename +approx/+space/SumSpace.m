classdef SumSpace < handle
    % SUMSPACE Sum of two spectral function spaces.
    %
    %   SumSpace represents the sum of exactly two spectral spaces:
    %
    %   \f[
    %       f(x) = \sum_{j=1}^{n_{lhs}} w_{lhs,j} c_{lhs,j} \phi_{lhs,j}(x)
    %            + \sum_{j=1}^{n_{rhs}} w_{rhs,j} c_{rhs,j} \phi_{rhs,j}(x)
    %   \f]
    %
    %   where \f$\phi_{lhs,j}(x)\f$ and \f$\phi_{rhs,j}(x)\f$ are basis
    %   functions from the left-hand side and right-hand side spaces,
    %   \f$c_{lhs,j}\f$ and \f$c_{rhs,j}\f$ are coefficients, and
    %   \f$w_{lhs,j}\f$ and \f$w_{rhs,j}\f$ are per-DOF weights.
    %
    % See also:
    %   approx.space.SpectralSpace

    properties
        Lhs % Left-hand side space
        Rhs % Right-hand side space
        LhsWeight % Weight vector for left-hand side space
        RhsWeight % Weight vector for right-hand side space
    end

    properties (Dependent)
        NDofs % Total number of degrees of freedom (integer)
        NDims % Number of spatial dimensions (integer)
        NLhsDofs % Number of DOFs in left-hand side space
        NRhsDofs % Number of DOFs in right-hand side space
        Weights % Combined weight vector for both spaces
    end

    methods
        function obj = SumSpace(lhs, rhs, options)
            % SUMSPACE Constructor for SumSpace.
            %
            %   obj = SumSpace(lhs, rhs) creates a sum of two spectral
            %   spaces. Both spaces must have the same spatial dimension.
            %
            %   obj = SumSpace(lhs, rhs, options) creates a sum with
            %   specified per-DOF weights. @a options is a struct with
            %   fields lhsWeight and rhsWeight.

            arguments
                lhs {mustBeA(lhs, {'numeric', 'approx.space.SpectralSpace'})}
                rhs {mustBeA(rhs, {'numeric', 'approx.space.SpectralSpace'})}
                options.lhsWeight {mustBeNumeric} = []
                options.rhsWeight {mustBeNumeric} = []
            end
            
            obj.Lhs = lhs;
            obj.Rhs = rhs;
            obj.setLhsWeight(options.lhsWeight);
            obj.setRhsWeight(options.rhsWeight);
        end

        function obj = setLhsWeight(obj, w)
            % SETLHSWEIGHT Set weight vector for left-hand side space.

            arguments
                obj approx.space.SumSpace
                w {mustBeNumeric}
            end

            core.except.assert(isempty(w) || length(w) == obj.Lhs.NDofs, ...
                'InvalidInput', ...
                'Weight vector length must be empty or match left-hand side DOFs.');

            obj.LhsWeight = w(:);
        end

        function obj = setRhsWeight(obj, w)
            % SETRHSWEIGHT Set weight vector for right-hand side space.

            arguments
                obj approx.space.SumSpace
                w {mustBeNumeric}
            end

            core.except.assert(isempty(w) || length(w) == obj.Rhs.NDofs, ...
                'InvalidInput', ...
                'Weight vector length must be empty or match right-hand side DOFs.');

            obj.RhsWeight = w(:);
        end

        function n = get.NDofs(obj)
            % GET.NDOFS Get the total number of degrees of freedom.

            n = obj.NLhsDofs + obj.NRhsDofs;
        end

        function n = get.NDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            if ~isempty(obj.Lhs)
                n = obj.Lhs.NDims;
                return;
            end

            if ~isempty(obj.Rhs)
                n = obj.Rhs.NDims;
                return;
            end

            n = 0;
        end

        function n = get.NLhsDofs(obj)
            % GET.NLHSDOFS Get the number of DOFs in left-hand side space.

            if isempty(obj.Lhs)
                n = 0;
            else
                n = obj.Lhs.NDofs;
            end
        end

        function n = get.NRhsDofs(obj)
            % GET.NRHSDOFS Get the number of DOFs in right-hand side space.

            if isempty(obj.Rhs)
                n = 0;
            else
                n = obj.Rhs.NDofs;
            end
        end

        function w = get.Weights(obj)
            % GET.WEIGHTS Get combined weight vector for both spaces.

            w = [obj.LhsWeight; obj.RhsWeight];
        end

        function Y = eval(obj, X, CLhs, CRhs)
            % EVAL Evaluate functions in the sum space.
            %
            %   Y = eval(obj, X, CLhs, CRhs) evaluates the sum of functions
            %   with coefficient matrix @a CLhs for the left-hand side space
            %   and coefficient matrix @a CRhs for the right-hand side space
            %   at evaluation points @a X.

            arguments
                obj approx.space.SumSpace
                X{mustBeNumeric}
                CLhs{mustBeNumeric}
                CRhs{mustBeNumeric}
            end

            YLhs = obj.evalLhs(X, CLhs);
            YRhs = obj.evalRhs(X, CRhs);
            Y = YLhs + YRhs;
        end

        function YLhs = evalLhs(obj, X, CLhs)
            % EVALLHS Evaluate functions in the left-hand side space only.
            %
            %   YLhs = evalLhs(obj, X, CLhs) evaluates functions with
            %   coefficient matrix @a CLhs for the left-hand side space
            %   at evaluation points @a X.

            arguments
                obj approx.space.SumSpace
                X{mustBeNumeric}
                CLhs{mustBeNumeric}
            end

            core.except.assert(~isempty(obj.Lhs), ...
                'InvalidOperation', ...
                'Left-hand side space is not defined.');

            DLhs = reshape(CLhs, size(CLhs, 1), []);
            if ~isempty(obj.LhsWeight)
                DLhs = obj.LhsWeight .* DLhs;
            end
            YLhs = obj.Lhs.eval(X, DLhs);
        end

        function YRhs = evalRhs(obj, X, CRhs)
            % EVALRHS Evaluate functions in the right-hand side space only.
            %
            %   YRhs = evalRhs(obj, X, CRhs) evaluates functions with
            %   coefficient matrix @a CRhs for the right-hand side space
            %   at evaluation points @a X.

            arguments
                obj approx.space.SumSpace
                X{mustBeNumeric}
                CRhs{mustBeNumeric}
            end

            core.except.assert(~isempty(obj.Rhs), ...
                'InvalidOperation', ...
                'Right-hand side space is not defined.');

            DRhs = reshape(CRhs, size(CRhs, 1), []);
            if ~isempty(obj.RhsWeight)
                DRhs = obj.RhsWeight .* DRhs;
            end
            YRhs = obj.Rhs.eval(X, DRhs);
        end

        function R = invLhs(obj, Y, X, CRhs)
            % INVLHS Extract left-hand side component by removing right-hand side.
            %
            %   R = invLhs(obj, Y, X, CRhs) extracts the left-hand side
            %   component by removing the right-hand side contribution from @a Y
            %   using coefficient matrix @a CRhs for the right-hand side space.

            arguments
                obj approx.space.SumSpace
                Y{mustBeNumeric}
                X{mustBeNumeric}
                CRhs{mustBeNumeric}
            end

            core.except.assert(~isempty(obj.Rhs), ...
                'InvalidOperation', ...
                'Right-hand side space is not defined.');

            YRhs = obj.Rhs.eval(X, CRhs);
            R = (Y - YRhs);
            if ~isempty(obj.LhsWeight)
                R = R ./ obj.LhsWeight(:).';
            end
        end

        function R = invRhs(obj, Y, X, CLhs)
            % INVRHS Extract right-hand side component by removing left-hand side.
            %
            %   R = invRhs(obj, Y, X, CLhs) extracts the right-hand side
            %   component by removing the left-hand side contribution from @a Y
            %   using coefficient matrix @a CLhs for the left-hand side space.

            arguments
                obj approx.space.SumSpace
                Y{mustBeNumeric}
                X{mustBeNumeric}
                CLhs{mustBeNumeric}
            end

            core.except.assert(~isempty(obj.Lhs), ...
                'InvalidOperation', ...
                'Left-hand side space is not defined.');

            YLhs = obj.Lhs.eval(X, CLhs);
            R = (Y - YLhs);
            if ~isempty(obj.RhsWeight)
                R = R ./ obj.RhsWeight(:).';
            end
        end

        function [CLhs, CRhs] = split(obj, C)
            % SPLIT Split coefficients into left and right components.
            %
            %   [CLhs, CRhs] = split(obj, C) splits the coefficient matrix
            %   @a C into left-hand side coefficients @a CLhs and
            %   right-hand side coefficients @a CRhs based on the DOFs
            %   of each space.

            arguments
                obj approx.space.SumSpace
                C{mustBeNumeric}
            end

            C = reshape(C, [], size(C, ndims(C)));
            nLhs = obj.NLhsDofs;

            CLhs = C(1:nLhs, :);
            CRhs = C(nLhs+1:end, :);
        end

        function C = combine(obj, CLhs, CRhs)
            % COMBINE Combine left and right coefficients.
            %
            %   C = combine(obj, CLhs, CRhs) combines the coefficient
            %   matrices @a CLhs and @a CRhs into a single coefficient
            %   matrix @a C.

            arguments
                obj approx.space.SumSpace
                CLhs{mustBeNumeric}
                CRhs{mustBeNumeric}
            end

            CLhs = reshape(CLhs, [], size(CLhs, ndims(CLhs)));
            CRhs = reshape(CRhs, [], size(CRhs, ndims(CRhs)));

            C = [CLhs; CRhs];
        end
    end
end