function positions = computeFigurePositions(count, screen)
%COMPUTEFIGUREPOSITIONS Compute bounded pixel positions without root mutations.
arguments
    count (1, 1) double {mustBePositive, mustBeInteger, mustBeFinite}
    screen (1, 4) double {mustBeFinite}
end
if any(screen(3:4) <= 0)
    error('mdemat:tools:visual:InvalidScreen', 'Screen dimensions must be positive.');
end
columns = ceil(sqrt(count)); rows = ceil(count / columns);
if columns > screen(3) || rows > screen(4)
    error('mdemat:tools:visual:FigureBudgetExceeded', ...
        'The screen cannot hold the requested number of positive-sized windows.');
end
margin = max(0, min([20, (screen(3) - columns) / (columns + 1), ...
    (screen(4) - rows) / (rows + 1)]));
width = max(1, min(560, (screen(3) - (columns + 1) * margin) / columns));
height = max(1, min(420, (screen(4) - (rows + 1) * margin) / rows));
positions = zeros(count, 4);
for k = 1:count
    column = mod(k - 1, columns); row = floor((k - 1) / columns);
    x = screen(1) + margin + column * (width + margin);
    y = screen(2) + screen(4) - margin - (row + 1) * height - row * margin;
    positions(k, :) = [x, max(screen(2), y), width, height];
end
end
