%==========================================================================
% FileName: ex01_periodic
% Author: Yi Cai
% Date: 2025-07-13
% Description:
%   Solving advection equation with periodic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 2; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ampl = [0, 1]; % Wave amplitude [sine, cosine]
config.ic = u0(config.omega, config.ampl); % Initial condition
config.bc = []; % Boundary condition
config.uExact = uExact(config.advection, config.omega, config.ampl); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.records = {'Integration', 'Visualization', 'Analysis'};
config.verbose = 1; % Verbose flag
config.TID = upper(config.tDiscName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('PA%d', config.nDims); % Experiment ID
config.SID = strjoin({config.TID, config.XID}, '-'); % Scheme ID

%% SIMULATION SETUP

timer = core.chrono.Timer(config.records, config.verbose);

visualizer = profilers.Visualizer(config.nDims);
visualizer.setDensity(ones(1, config.nDims));
visualizer.setTimeline(config.tFinal, 1);
visualizer.setComponents(struct('U', 1));
visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
visualizer.addExact('U', config.uExact);
visualizer.addDataset('REF', 1);
visualizer.addDataset(config.XID, 2);

switch config.nDims
    case 1
        visualizer.addPlotter('profile', '1d');
    case 2
        visualizer.addPlotter('profile', '2d');
        visualizer.addPlotter('slice', 'slice1d');
    case 3
        visualizer.addPlotter('profile', '3d');
        visualizer.addPlotter('slice', 'slice1d');
end

analyzer = profilers.Analyzer(5);
analyzer.setReductions({'L1', 'L2', 'Lx'});
analyzer.setDensity(repmat(4, 1, config.nDims));
analyzer.setComponents(struct('U', 1));
analyzer.addExact('U', config.uExact);
analyzer.initialize();

if contains(config.XID, 'sl')
    cls = 'SldgScheme';
else
    cls = 'UwdgScheme';
end

scheme = physics.advection.(cls)(config);
scheme.setTimeDiscretization(config.tDiscName, config.tFinal);

%% SIMULATION EXECUTION

xElement = approx.element.L2Element.orthotope(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);

if ~analyzer.hasExact
    %< Mesh initialization
    xMesh = approx.mesh.Grid.uniform(config.xBBox, config.nx/2);
    for i = 1:analyzer.nLevels + 1
        %< Reset timer
        timer.reset();

        %< Create state
        xDisc = approx.space.MeshSpace(xElement, xMesh);
        state = physics.State(xDisc);

        %< Initialization
        state = scheme.initialize(state);
        if visualizer.isEnabled
            visualizer.initialize();
            timer.start(2, '[M] Visualization.');
            visualizer.plot(state.xDisc, state.dofs, scheme.tDisc);
            timer.stop(2);
        end

        %< Simulation
        while ~scheme.isExhausted
            state = scheme.preStep(state);

            timer.start(1, scheme.echo(timer.records(1).duration));
            state = scheme.step(state);
            timer.stop(1);

            if visualizer.isEnabled
                visualizer.initialize();
                timer.start(2, '[M] Visualization.');
                visualizer.plot(state.xDisc, state.dofs, scheme.tDisc);
                timer.stop(2);
            end

            state = scheme.postStep(state);
        end

        %< Finalization
        state = scheme.finalize(state);

        %< Error evaluation
        if i > 1 && analyzer.isEnabled
            timer.start(3, '[M] Error Evaluation.');
            analyzer.addLevel(xMesh);
            analyzer.addRelative(state.xDisc, state.dofs, coarse.xDisc, coarse.dofs);
            timer.stop(3);
        end

        %< Save coarse state
        coarse = state;

        %< Mesh refinement
        xMesh = xMesh.refine(1);
    end
else
    %< Mesh initialization
    xMesh = approx.mesh.Grid.uniform(config.xBBox, config.nx);
    for i = 1:analyzer.nLevels
        %< Reset timer
        timer.reset();

        %< Create state
        xDisc = approx.space.MeshSpace(xElement, xMesh);
        state = physics.State(xDisc);

        %< Initialization
        state = scheme.initialize(state);

        if visualizer.isEnabled
            visualizer.initialize();
            timer.start(2, '[M] Visualization.');
            visualizer.plot(state.xDisc, state.dofs, scheme.tDisc);
            timer.stop(2);
        end

        %< Simulation
        while ~scheme.isExhausted
            state = scheme.preStep(state);

            timer.start(1, scheme.echo(timer.records(1).duration));
            state = scheme.step(state);
            timer.stop(1);

            if visualizer.isEnabled
                timer.start(2, '[M] Visualization.');
                visualizer.plot(state.xDisc, state.dofs, scheme.tDisc);
                timer.stop(2);
            end

            state = scheme.postStep(state);
        end

        %< Finalization
        state = scheme.finalize(state);

        %< Error evaluation
        if analyzer.isEnabled
            timer.start(3, '[M] Error Evaluation.');
            analyzer.addLevel(xMesh);
            analyzer.addAbsolute(state.xDisc, state.dofs, config.tFinal);
            timer.stop(3);
        end

        %< Mesh refinement
        xMesh = xMesh.refine(1);
    end
end

if config.verbose > 0
    timer.report();
end

results = analyzer.analyze(config.verbose);
% analyzer.draw(repmat(config.xBasisOrder, analyzer.nReductions));

%% HELPER FUNCTIONS
function h = u0(omega, ampl)
h = @(x) u0Impl(x, omega, ampl);

    function f = u0Impl(x, omega, ampl)
    z = omega(:).' * x;
    f = ampl(1) * cos(z) + ampl(2) * sin(z);
    end
end

function h = uExact(c, omega, ampl)
h = @(x, t) uExactImpl(x, t, c, omega, ampl);
    function f = uExactImpl(x, t, c, omega, ampl)
        z = omega(:).' * (x - c(:) * t);
        f = ampl(1) * cos(z) + ampl(2) * sin(z);
    end
end