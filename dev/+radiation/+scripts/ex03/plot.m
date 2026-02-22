clc;
clear;

M = 1;
N = 8;
EXPERIMENT_ID = 'DSID1R2';
XDISC_ID = 'UWLDG';
DATA_DIR = fullfile(fileparts(mfilename('fullpath')));
DEFAULT_MARKERS = {'o', 's', '^', 'v', 'd', 'x', '+'};
DEFAULT_COLORS = {'#000000', '#004488', '#BB5566', '#DDAA33'};

dataset = struct();
style = struct();

file = sprintf('%s-REF.mat', EXPERIMENT_ID);
path = fullfile(DATA_DIR, file);
if exist(path, 'file')
    data = load(path).REF;
    dataset.REF = data;
    style.REF = { ...
        'Color', DEFAULT_COLORS{1}, ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        };
end

P = { ...
    'BE', 1, M, N; ...
    'SDIRK2', 2, M, N; ...
    'SDIRK3', 3, M, N; ...
    };
n = 0;
for i = 1:size(P, 1)
    schemeId = sprintf('%s-%s%d-P%dS%d', P{i, 1}, XDISC_ID, ...
        P{i, 2}, P{i, 3}, P{i, 4});
    schemeName = sprintf('%s%d', XDISC_ID, P{i, 2});
    file = sprintf('%s.mat', strjoin({EXPERIMENT_ID, schemeId}, '-'));
    path = fullfile(DATA_DIR, file);
    if exist(path, 'file')
        n = n + 1;
        data = load(path).DG;
        dataset.(schemeName) = data;
        style.(schemeName) = { ...
            'Color', DEFAULT_COLORS{1+n}, ...
            'Marker', DEFAULT_MARKERS{n}, ...
            'LineStyle', 'none', ...
            'LineWidth', 1, ...
            };
    end
end


figure(1);

subdata = prepareData(dataset, 'U1', 1);


% scheme.visualizer.reset();
% scheme.visualizer.render([], dataset, style, []);
% ax = gca;
% set(ax, 'box', 'off');
% ax.XAxis.TickValues = 0:0.1:1;
% ax.XAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues = 0:0.02:1;
% ax.YAxis.TickValues = -0.2:0.2:1.4;
% ax.YAxis.MinorTick = 'on';
% ax.YAxis.MinorTickValues = -0.2:0.04:1.4;
% xlim([0, 1]);
% ylim([-0.2, 1.4]);
% xl = xlabel("$x$", "Interpreter", "latex");
% yl = ylabel("$\rho$", "Interpreter", "latex");
% xl.FontSize = 16;
% yl.FontSize = 16;