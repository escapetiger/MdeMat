classdef FourierBasisCompiler < core.symbolic.Compiler
    % FOURIERBASISCOMPILER Compiler for Fourier basis functions.
    %
    %   FourierBasisCompiler specializes in compiling basis matrices for
    %   real-valued Fourier basis functions. It creates and manages arrays
    %   of Fourier basis functions including constant, cosine, and sine
    %   terms on specified intervals.
    %
    % See Also:
    %   core.symbolic.Compiler, core.symbolic.FourierBasisFunction
    
    methods
        function obj = setFunctions(obj, maxIndex, bbox, basisType)
            % SETFUNCTIONS Sets Fourier basis functions.
            %
            %   obj = setFunctions(obj, maxIndex) creates an array of
            %   Fourier basis functions from index 1 to maxIndex on the
            %   default interval [0, 2π].
            %   
            %   obj = setFunctions(obj, maxIndex, bbox) creates Fourier
            %   basis functions on the specified interval.
            %   
            %   obj = setFunctions(obj, maxIndex, bbox, basisType) creates
            %   Fourier basis functions with explicit basis type
            %   specification.
            %
            % Inputs:
            %   obj - The FourierBasisCompiler object
            %   maxIndex - Maximum index for basis functions (positive integer)
            %   bbox - Bounding box [lower, upper] (optional, default:[0,2*pi])
            %   basisType - Basis type (optional, default:'fourier')
            %
            % Outputs:
            %   obj - The FourierBasisCompiler object
            
            core.except.assert(maxIndex >= 1 && isscalar(maxIndex) && ...
                maxIndex == floor(maxIndex), 'InvalidInput', ...
                'Maximum index must be a positive integer.');
            
            if nargin < 3 || isempty(bbox)
                bbox = [0, 2*pi];
            end
            
            core.except.assert(isvector(bbox) && length(bbox) == 2, ...
                'InvalidInput', 'Bounding box must be [lower, upper].');
            
            core.except.assert(bbox(1) < bbox(2), ...
                'InvalidInput', 'Lower bound must be less than upper bound.');
            
            if nargin < 4
                basisType = 'fourier';
            end
            
            supportedTypes = {'fourier'};
            core.except.assert(ismember(lower(basisType), supportedTypes), ...
                'InvalidInput', 'Supported basis types: %s', strjoin(supportedTypes, ', '));
            
            x = obj.variables;
            
            obj.funcs = arrayfun(@(i) core.symbolic.FourierBasisFunction(i, bbox, x), ...
                1:maxIndex);
            
            fprintf('Set %d %s basis functions, index 1 to %d on [%.6g, %.6g]\n', ...
                maxIndex, basisType, maxIndex, bbox(1), bbox(2));
        end
    end
end