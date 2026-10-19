classdef Metric < handle
    % METRIC Base class for error and distance metrics.
    %
    %   Metric defines the interface for computing distances and errors
    %   between abstract mathematical objects such as functions, solutions,
    %   or data sets. This class serves as the foundation for all specific
    %   metric implementations used in numerical analysis and error
    %   computation.
    %
    %   All concrete metric classes must inherit from this base class and
    %   implement the abstract evaluate method to define how errors are
    %   computed for their specific use case. The class supports various
    %   error reduction types commonly used in numerical analysis.
    
    properties (Constant)
        VALID_REDUCTION = {'L1', 'L2', 'Lx', 'none'} % Valid error reduction types
    end

    properties
        reduction % Cell array of reduction types to compute
    end

    methods
        function obj = Metric(reduction)
            % METRIC Constructor for Metric base class.
            %
            %   obj = Metric(reduction) creates a metric with the specified
            %   error reduction types. The reduction parameter determines
            %   which error norms will be computed when evaluate is called.
            %
            % Inputs:
            %   reduction - Cell array of reduction types:
            %
            % Notes:
            %   reduction = 'L1': L_1 norm (sum of absolute values)
            %   reduction = 'L2': L_2 norm (Euclidean norm)
            %   reduction = 'Lx': L_inf norm (maximum absolute value)
            %   reduction = 'none': No reduction (return full error array)
            %
            % Outputs:
            %   obj - Constructed Metric object

            %< Validate reduction types
            if ischar(reduction)
                reduction = {reduction};
            end
            
            for i = 1:length(reduction)
                core.except.assert(ismember(reduction{i}, obj.VALID_REDUCTION), ...
                    'InvalidInput', 'Invalid reduction type: %s', reduction{i});
            end

            obj.reduction = reduction;
        end
    end

    methods (Abstract)
        % EVALUATE Evaluate error or distance between objects.
        error = evaluate(obj, varargin)
    end
end