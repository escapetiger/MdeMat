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
    
    properties
        Reductions % Cell array of reduction types to compute
    end

    methods
        function obj = Metric(options)
            % METRIC Construct an instance of Metric.
            %
            %   obj = Metric() creates a metric with empty reductions list.
            %
            %   obj = Metric(reductions=reductions) creates a metric with
            %   the specified @a reductions cell array of error reduction
            %   types.
            %
            % See also:
            %   physics.analysis.FiniteElementAbsoluteMetric
            
            arguments
                options.reductions cell = {}
            end
            
            obj.Reductions = options.reductions;
        end
        
        function obj = addReduction(obj, reduction)
            % ADDREDUCTION Add error reduction type to the metric.
            %
            %   obj = addReduction(obj, reduction) adds the specified
            %   error reduction type @a reduction to the metric's list of
            %   reductions to compute during error evaluation.
            
            arguments
                obj physics.analysis.Metric
                reduction {mustBeTextScalar}
            end
            
            obj.Reductions{end+1} = reduction;
        end
    end
end 