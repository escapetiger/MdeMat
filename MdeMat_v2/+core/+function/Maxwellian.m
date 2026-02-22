classdef Maxwellian < core.function.Function
    % MAXWELLIAN Shifted Maxwell-Boltzmann distribution function.
    %
    %   Maxwellian represents a shifted Maxwell-Boltzmann distribution:
    %
    %   \f[
    %     f(x;y,T) = (1/(2*pi*T)^{d/2}) * exp(-|x-y|^2/(2*T))
    %   \f]
    %
    %   where \f$T\f$ is the temperature, \f$y\f$ is the center point, and
    %   \f$d\f$ is the dimension. This represents a Gaussian distribution
    %   centered at point y with covariance matrix T*I.
    %
    % Examples:
    %   % 1D Maxwellian centered at y=2 with T T=1
    %   f = MaxwellianFunction(1, 1, 2);  % R^1 -> R^1
    %
    %   % 3D velocity distribution centered at drift velocity
    %   y = [1; 0; -0.5];  % Drift velocity
    %   f = MaxwellianFunction(3, 0.5, y);  % R^3 -> R^1
    %
    %   % Evaluate at points
    %   X = randn(3, 1000);
    %   Y = f.evaluate(X);
    %
    % See Also:
    %   core.function.Function

    properties
        temperature % Temperature
        center % Center point (nDims x 1 vector)
    end

    methods
        function obj = Maxwellian(nDims, temperature, center)
            % MAXWELLIAN Constructor for the Maxwellian
            % class.
            %
            %   obj = Maxwellian(nDims, temperature, center) creates a
            %   Maxwellian in nDims dimensions with the specified
            %   temperature and center point.
            %
            % Inputs:
            %   nDims - Number of input dimensions (positive integer)
            %   temperature - Temperature parameter
            %   center - Center point (nDims x 1 vector or scalar)
            %
            % Outputs:
            %   obj - The created Maxwellian object

            if nargin < 3 || isempty(center)
                center = zeros(nDims, 1);
            end

            obj@core.function.Function(nDims, 1);
            obj.temperature = temperature;
            obj.center = center(:);
        end
    end

    methods (Access = protected)
        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of Maxwellian evaluation.
            %
            %   Y = evaluateImpl(obj, X) evaluates the Maxwellian
            %   distribution at the specified points.
            %
            % Inputs:
            %   obj - The MaxwellianFunction object
            %   X - Coordinates (nDims x nPoints)
            %
            % Outputs:
            %   Y - Function values (1 x nPoints)

            d = obj.nDims;
            T = obj.temperature;
            C = obj.center;
            Y = exp(-sum((X - C).^2, 1) ./ (2 * T)) ./ (2 * pi * T).^(d/2);
        end

        function dY = derivativeImpl(obj, X, r)
            % DERIVATIVEIMPL Implementation of Maxwellian derivative.
            %
            %   dY = derivativeImpl(obj, X, r) computes derivatives of the
            %   Maxwellian function.
            %
            %   For first derivatives: 
            %   \f[
            %     \frac{\partial{f}}{\partial{x_i}} = -(x_i-y_i)/T * f(x)
            %   \f]
            %   For higher derivatives: not implemented
            %
            % Inputs:
            %   obj - The Maxwellian object
            %   X - Coordinates (nDims x nPoints)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (1 x nPoints)

            T = obj.temperature;
            if sum(r) == 1
                Y = obj.evaluateImpl(X);
                dX = X - obj.center;
                dY = -dX(find(r > 0, 1), :) / T .* Y;
            else
                core.except.assert(0, 'InvalidInput', ...
                    'High-order derivatives are not supported.');
            end
        end
    end
end