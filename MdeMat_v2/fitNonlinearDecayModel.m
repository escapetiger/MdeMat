function fitResult = fitNonlinearDecayModel()
% fitNonlinearDecayModel
% Fit data using the nonlinear decay model y = a - b*(x-c)^(-d).
%
% Syntax
%   fitResult = fitNonlinearDecayModel()
%
% Description
% This function fits experimental data to a nonlinear decay model of the
% form y = a - b*(x-c)^(-d), where a, b, c, and d are parameters to be estimated.
% The parameter 'a' is constrained to be in the range [9, 12].
% The function uses nonlinear least squares fitting and displays results
% with visualization.
%
% Output Arguments
%   fitResult - structure containing fitted parameters and statistics
%     * parameters - fitted values of a, b, c, and d
%     * rsquared - coefficient of determination
%     * rmse - root mean square error
%     * residuals - fitting residuals
%
% Example
%   result = fitNonlinearDecayModel()
%
% See also
% fitnlm, nlinfit, lsqcurvefit

%% Initialize data
% X-axis values
xData = [0, 1, 2, 3, 4, 5, 8, 10, 12, 15, 18, 20, 25, 30, 35, 40, 50, 55, 60, 70, 72, 74, 80]';

% Y-axis values  
yData = [0, 6.2097, 6.7417, 7.2635, 7.5920, 7.8747, 8.4854, 8.7461, 8.9925, ...
         9.2179, 9.3467, 9.4431, 9.6278, 9.7492, 9.8206, 9.8540, 9.9247, ...
         9.9654, 9.9722, 10.0169, 10.0064, 10.0092, 10.0231]';

% Validate input data
validateInputData(xData, yData);

%% Define model function
% Model: y = a - b*(x-c)^(-d)
modelFunction = @(parameters, x) calculateModelValues(parameters, x);

%% Set initial parameter estimates
% Based on the data pattern, provide reasonable initial guesses
initialA = 10.0;   % Initial guess for asymptote parameter a (within [9,12])
initialB = 10.0;   % Initial guess for parameter b
initialC = -1.0;   % Initial guess for shift parameter c
initialD = 0.5;    % Initial guess for power parameter d (positive for negative exponent)
initialParameters = [initialA, initialB, initialC, initialD];

%% Perform nonlinear fitting
try
    % Use lsqcurvefit for robust nonlinear fitting
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 1000, ...
                          'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8);
    
    % Set parameter bounds to ensure physical meaningfulness and constraints
    lowerBounds = [9, 0, -Inf, 0];    % a ∈ [9,12], b > 0, c unconstrained, d > 0
    upperBounds = [12, Inf, Inf, Inf];   % Upper bound for a is 12
    
    [fittedParameters, residualNorm] = lsqcurvefit(modelFunction, ...
        initialParameters, xData, yData, lowerBounds, upperBounds, options);
    
catch fittingError
    error('fitNonlinearDecayModel:fittingFailed', ...
        'Nonlinear fitting failed: %s', fittingError.message);
end

%% Calculate fitted values and statistics
fittedY = modelFunction(fittedParameters, xData);
residuals = yData - fittedY;

% Calculate goodness-of-fit statistics
rSquared = calculateRSquared(yData, fittedY);
rmse = sqrt(mean(residuals.^2));

%% Display results
displayFittingResults(fittedParameters, rSquared, rmse);

%% Create visualization
createFittingPlot(xData, yData, fittedY, fittedParameters, rSquared);

%% Prepare output structure
fitResult = struct();
fitResult.parameters = struct('a', fittedParameters(1), 'b', fittedParameters(2), ...
                             'c', fittedParameters(3), 'd', fittedParameters(4));
fitResult.rsquared = rSquared;
fitResult.rmse = rmse;
fitResult.residuals = residuals;
fitResult.fittedValues = fittedY;

end

%% Helper Functions

function yValues = calculateModelValues(parameters, x)
% Calculate model values y = a - b*(x-c)^(-d) with proper handling of edge cases
%
% Input Arguments
%   parameters - array containing [a, b, c, d]
%   x - input x values
%
% Output Arguments
%   yValues - calculated y values

a = parameters(1);
b = parameters(2);
c = parameters(3);
d = parameters(4);

% Calculate (x-c)
xShifted = x - c;

% For negative exponents, we need to handle division by zero and negative bases
% Add small epsilon to avoid division by zero
epsilon = 1e-10;

% Ensure we don't have zero or negative bases for negative exponents
xShiftedSafe = sign(xShifted) .* max(abs(xShifted), epsilon);

% Calculate y = a - b*(x-c)^(-d) = a - b / (x-c)^d
yValues = a - b ./ (xShiftedSafe.^d);

% Replace any NaN or Inf values with a large residual value to penalize invalid solutions
invalidIndices = ~isfinite(yValues);
if any(invalidIndices)
    yValues(invalidIndices) = 1e6;  % Large value to penalize invalid solutions
end

end

function validateInputData(xData, yData)
% Validate input data for fitting
%
% Input Arguments
%   xData - x-axis data values
%   yData - y-axis data values

if length(xData) ~= length(yData)
    error('fitNonlinearDecayModel:dataMismatch', ...
        'X and Y data must have the same length.');
end

if any(isnan(xData)) || any(isnan(yData))
    error('fitNonlinearDecayModel:nanValues', ...
        'Data cannot contain NaN values.');
end

end

function rSquared = calculateRSquared(observedY, predictedY)
% Calculate coefficient of determination (R-squared)
%
% Input Arguments
%   observedY - observed y values
%   predictedY - predicted y values from model
%
% Output Arguments
%   rSquared - coefficient of determination

sumSquaresTotal = sum((observedY - mean(observedY)).^2);
sumSquaresResidual = sum((observedY - predictedY).^2);
rSquared = 1 - (sumSquaresResidual / sumSquaresTotal);

end

function displayFittingResults(fittedParameters, rSquared, rmse)
% Display fitting results in formatted output
%
% Input Arguments
%   fittedParameters - array containing fitted a, b, c, and d values
%   rSquared - coefficient of determination
%   rmse - root mean square error

fprintf('\n=== Nonlinear Fitting Results ===\n');
fprintf('Model: y = a - b*(x-c)^(-d)\n\n');
fprintf('Fitted Parameters:\n');
fprintf('  a = %.4f (constrained to [9, 12])\n', fittedParameters(1));
fprintf('  b = %.4f\n', fittedParameters(2));
fprintf('  c = %.4f\n', fittedParameters(3));
fprintf('  d = %.4f\n', fittedParameters(4));
fprintf('\nModel equation with fitted parameters:\n');
fprintf('  y = %.4f - %.4f*(x-%.4f)^(-%.4f)\n', ...
    fittedParameters(1), fittedParameters(2), fittedParameters(3), fittedParameters(4));
fprintf('\nGoodness of Fit:\n');
fprintf('  R-squared = %.4f\n', rSquared);
fprintf('  RMSE = %.4f\n', rmse);
fprintf('================================\n\n');

end

function createFittingPlot(xData, yData, fittedY, fittedParameters, rSquared)
% Create visualization of original data and fitted curve
%
% Input Arguments
%   xData - x-axis data values
%   yData - observed y values
%   fittedY - fitted y values from model
%   fittedParameters - fitted parameter values
%   rSquared - coefficient of determination

% Create new figure with appropriate size
figureHandle = figure('Position', [100, 100, 900, 600]);

% Plot original data points
scatter(xData, yData, 100, 'bo', 'filled', 'DisplayName', 'Original Data');
hold on;

% Generate smooth curve for fitted model
% Be careful with the range to avoid numerical issues
xMin = max(min(xData), fittedParameters(3) + 0.1);  % Avoid singularity at x = c
xMax = max(xData);
xSmooth = linspace(xMin, xMax, 300);
ySmooth = calculateModelValues(fittedParameters, xSmooth);

% Only plot valid points
validIndices = isfinite(ySmooth) & (ySmooth > -100) & (ySmooth < 100);
plot(xSmooth(validIndices), ySmooth(validIndices), 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');

% Plot fitted points for comparison
plot(xData, fittedY, 'go', 'MarkerSize', 6, 'DisplayName', 'Fitted Points');

% Enhance plot appearance
xlabel('X Values', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Values', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Nonlinear Fitting: y = a - b*(x-c)^{-d}\na = %.4f, b = %.4f, c = %.4f, d = %.4f, R² = %.4f', ...
    fittedParameters(1), fittedParameters(2), fittedParameters(3), ...
    fittedParameters(4), rSquared), 'FontSize', 12, 'FontWeight', 'bold');

% Add grid and legend
grid on;
legend('Location', 'southeast', 'FontSize', 11);

% Set axis properties
set(gca, 'FontSize', 11);
xlim([min(xData) - 2, max(xData) + 2]);
ylim([min(yData) - 0.5, max(yData) + 0.5]);

hold off;

end