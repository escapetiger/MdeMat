classdef FunctionMetric < profilers.metrics.Metric
    % FUNCTIONMETRIC Error computation for function-based solutions.
    %
    %   FunctionMetric implements error computation between discrete
    %   function values evaluated at specified integration nodes with
    %   associated quadrature weights. This class is designed for computing
    %   various error norms (L1, L2, Lx) between numerical and reference
    %   solutions represented as function values.
    %
    %   The metric uses numerical integration to compute weighted error
    %   norms, making it suitable for finite element error analysis where
    %   solutions are evaluated at quadrature points with appropriate
    %   integration weights.
    %
    % Examples:
    %   % Create function metric with L2 norm
    %   metric = FunctionMetric({'L2'}, gaussIntegrator);
    %   
    %   % Compute error between numerical and exact solutions
    %   numericalValues = [1.1, 2.05, 2.98]; % at integration points
    %   exactValues = [1.0, 2.0, 3.0];
    %   error = metric.evaluate(numericalValues, exactValues);
    %
    %   % Multiple error norms simultaneously
    %   metric = FunctionMetric({'L1', 'L2', 'Lx'}, integrator);
    %   errors = metric.evaluate(numerical, exact); % returns cell array
    %
    % See also:
    %   profilers.metrics.Metric, profilers.metrics.MeshSpaceAbsoluteMetric
    
    properties
        integrator % Numerical integrator providing nodes and weights
    end
    
    methods        
        function obj = FunctionMetric(reduction, integrator)
            % FUNCTIONMETRIC Constructor for FunctionMetric.
            %
            %   obj = FunctionMetric(reduction, integrator) creates a function
            %   metric that computes errors between function values using
            %   the specified integration rule for weighted error norms.
            %
            % Inputs:
            %   reduction - Cell array of error reduction types
            %   integrator - Integrator object providing nodes and weights
            %
            % Outputs:
            %   obj - Constructed FunctionMetric object

            obj@profilers.metrics.Metric(reduction);
            obj.integrator = integrator;
        end
        
        function E = evaluate(obj, U, V)
            % EVALUATE Compute error between function-based solutions.
            %
            %   E = evaluate(obj, U, V) computes errors between numerical
            %   solution U and reference solution V using the specified
            %   error reduction types. For L1 and L2 norms, integration
            %   weights are applied for proper numerical integration.
            %
            % Inputs:
            %   obj - The FunctionMetric object
            %   U - Numerical solution values at integration points
            %   V - Reference solution values (same points) or scalar
            %
            % Outputs:
            %   E - Cell array of computed errors for each reduction type
            
            %< Prepare data for weighted integration
            if ismember('L1', obj.reduction) || ismember('L2', obj.reduction)
                I = obj.integrator;
                w = I.weights;
                n = I.nPoints;
                U = reshape(U, n, []);
                if ~isscalar(V)
                    V = reshape(V, n, []);
                end
            end
            
            %< Compute absolute difference
            R = abs(U - V);
            
            %< Compute each requested error norm
            nReds = length(obj.reduction);
            E = cell(1, nReds);
            for i = 1:nReds
                red = obj.reduction{i};
                switch red
                    case 'L1'
                        e = sum(w * R, 'all');
                    case 'L2'
                        e = sqrt(sum(w * R.^2, 'all'));
                    case 'Lx'
                        e = max(R, [], 'all');
                    case 'none'
                        e = R;
                end
                E{i} = e(:).';
            end
        end
    end
end