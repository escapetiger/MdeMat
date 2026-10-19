classdef CompiledSeparableBasis < core.function.SeparableFunction
    % COMPILEDSEPARABLEBASIS Separable function composed of basis function
    % factors with optional sparse structure.
    %
    %   CompiledSeparableBasis represents a separable function where each
    %   factor is either a CompiledInterpolationBasis or a
    %   CompiledOrthogonalBasis. This provides a unified framework for
    %   creating tensor product basis functions with optional sparsity
    %   patterns for high-dimensional approximation.
    %
    %   The class supports both dense (full tensor product) and sparse
    %   (reduced set of multi-indices) evaluation patterns for both
    %   function values and derivatives, with automatic metadata loading
    %   from the compiled basis functions.
    %
    % Examples:
    %   % Dense separable basis with Legendre polynomials
    %   factors = {CompiledOrthogonalBasis(5, 'legendre', 'canonical'), ...
    %              CompiledOrthogonalBasis(4, 'legendre', 'canonical')};
    %   basis = CompiledSeparableBasis(factors);
    %   basis = basis.autoLoad();
    %   
    %   % Sparse basis with L1 constraint using Lagrange interpolation
    %   factors = {CompiledInterpolationBasis(6, 'lagrange', 'unit', 'gauss_legendre'), ...
    %              CompiledInterpolationBasis(5, 'lagrange', 'unit', 'gauss_legendre')};
    %   basis = CompiledSeparableBasis(factors, 'l1');
    %   basis = basis.autoLoad();
    %   
    %   % Mixed basis types with evaluation
    %   factors = {CompiledOrthogonalBasis(4, 'legendre', 'unit'), ...
    %              CompiledInterpolationBasis(3, 'hermite', 'canonical', 'gauss_lobatto')};
    %   basis = CompiledSeparableBasis(factors);
    %   basis = basis.autoLoad();
    %   y = basis.evaluate([0.5; 0.3]);
    %   dy = basis.differentiate([0.5; 0.3], [1, 0]);
    %
    % See Also:
    %   core.function.SeparableFunction,
    %   approx.basis.CompiledInterpolationBasis,
    %   approx.basis.CompiledOrthogonalBasis

    properties (Dependent)
        hasMetadata % True if all factors have loaded metadata
        metasource  % Cell array of available metadata sources for all factors
    end

    methods
        function obj = autoLoad(obj)
            % AUTOLOAD Automatically load metadata for all basis factors.
            %
            %   obj = autoLoad(obj) calls autoLoad() on each factor
            %   function to load appropriate compiled metadata. This
            %   enables efficient evaluation of basis functions and their
            %   derivatives by loading precompiled function handles for all
            %   factors.
            %
            % Inputs:
            %   obj - The CompiledSeparableBasis object
            %
            % Outputs:
            %   obj - Modified object with loaded metadata for all factors

            F = obj.factors;
            if isempty(F)
                return;
            end

            if iscell(F)
                for i = 1:length(F)
                    F{i} = F{i}.autoLoad();
                end
                obj.setFactors(F);
            else
                for i = 1:length(F)
                    F(i) = F(i).autoLoad();
                end
                obj.setFactors(F);
            end
        end

        function obj = load(obj, filenames)
            % LOAD Load metadata from specified files for all factors.
            %
            %   obj = load(obj, filenames) loads compiled metadata from the
            %   specified filenames for each factor function. The number of
            %   filenames must match the number of factors. This provides
            %   manual control over which metadata files are loaded.
            %
            % Inputs:
            %   obj - The CompiledSeparableBasis object
            %   filenames - Metadata files for each factor
            %
            % Outputs:
            %   obj - Modified object with loaded metadata for all factors

            F = obj.factors;
            if isempty(F)
                return;
            end

            core.except.assert(iscell(filenames), ...
                'InvalidInput', ...
                'Filenames must be a cell array.');

            core.except.assert(length(filenames) == length(F), ...
                'InvalidInput', ...
                'Number of filenames (%d) must match number of factors (%d).', ...
                length(filenames), length(F));

            if iscell(F)
                for i = 1:length(F)
                    F{i} = F{i}.load(filenames{i});
                end
                obj.setFactors(F);
            else
                for i = 1:length(F)
                    F(i) = F(i).load(filenames{i});
                end
                obj.setFactors(F);
            end
        end

        function TF = get.hasMetadata(obj)
            % GET.HASMETADATA Getter for hasMetadata dependent property.

            F = obj.factors;
            if isempty(F)
                TF = true;
                return;
            end
            
            if iscell(F)
                TF = all(cellfun(@(f) f.hasMetadata, F));
            else
                TF = all(arrayfun(@(f) f.hasMetadata, F));
            end
        end

        function sources = get.metasource(obj)
            % GET.METASOURCE Getter for metasource dependent property.

            F = obj.factors;
            if isempty(F)
                sources = {};
                return;
            end

            if iscell(F)
                sources = cellfun(@(f) f.metasource, F, 'Un', 0);
            else
                sources = arrayfun(@(f) f.metasource, F, 'Un', 0);
            end
        end
    end
end