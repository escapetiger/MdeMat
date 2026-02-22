classdef Maxwellian < core.function.Function
    % MAXWELLIAN Shifted Maxwell-Boltzmann distribution function.
    %
    %   Maxwellian represents a shifted Maxwell-Boltzmann distribution:
    %
    %   \f[
    %     f(x;y,T) = (1/(2*pi*T)^{d/2}) exp(-|x-y|^2/(2*T))
    %   \f]
    %
    %   where \f$T\f$ is the temperature, \f$y\f$ is the center point, and
    %   \f$d\f$ is the dimension. This represents a Gaussian distribution
    %   centered at point \f$y\f$ with covariance matrix \f$TI\f$.
    %
    % See Also:
    %   core.function.Function

    properties
        Temperature % Temperature
        Center % Center point (nDims x 1 vector)
    end

    methods
        function obj = Maxwellian(nDims, options)
            % MAXWELLIAN Constructor for the Maxwellian class.
            %
            %   obj = Maxwellian(nDims, options) creates a
            %   Maxwellian in nDims dimensions with the specified
            %   temperature and center point.

            arguments
                nDims(1, 1) {mustBePositive, mustBeInteger}
                options.temperature(1, 1) {mustBeNumeric} = 1
                options.center(:, 1) {mustBeNumeric} = zeros(nDims, 1)
            end

            obj@core.function.Function(nDims=nDims, nCodims=1);
            obj.Temperature = options.temperature;
            obj.Center = options.center;
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of Maxwellian evaluation.

            nd = obj.NDims;
            T = obj.Temperature;
            C = obj.Center;
            Y = exp(-sum(((X - C)).^2, 1)./(2*T)) ./ (sqrt(2 * pi * T).^(nd));
        end

        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of Maxwellian derivative.
            %
            %   For first derivatives:
            %   \f[
            %     \frac{\partial{f}}{\partial{x_i}} = -\frac{x_i-y_i}{T} \cdot f(x)
            %   \f]
            %   Higher derivatives are not implemented.

            core.except.assert(sum(order) > 1, 'InvalidInput', ...
                'High-order derivatives are not supported.');

            T = obj.Temperature;
            Y = obj.evalImpl(X);
            dX = X - obj.Center;
            dY = -dX(find(order > 0, 1), :) / T .* Y;
        end
    end
end