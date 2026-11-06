function style = plotStyle(options)
%PLOTSTYLE Create an independent value record of scientific plotting options.
arguments
    options.discreteColorScheme (1, 1) string = "bright"
    options.continuousColorScheme (1, 1) string = "jet"
    options.lineScheme cell = {'-'}
    options.scatterScheme cell = {'o', 's', '^', 'v', 'd', 'x', '+'}
    options.lineWidth (1, 1) double {mustBePositive, mustBeFinite} = 2
    options.fontSize (1, 1) double {mustBePositive, mustBeFinite} = 18
    options.axisLabelFontSize (1, 1) double {mustBePositive, mustBeFinite} = 28
    options.legendFontSize (1, 1) double {mustBePositive, mustBeFinite} = 18
    options.titleFontSize (1, 1) double {mustBePositive, mustBeFinite} = 28
    options.legendLocation (1, 1) string = "best"
    options.colorPadding (1, 1) double {mustBeNonnegative, mustBeFinite} = 0.05
    options.symmetricColorLimits (1, 1) logical = false
    options.logX (1, 1) logical = false
    options.logY (1, 1) logical = false
    options.complexDisplay (1, 1) string {mustBeMember(options.complexDisplay, ...
        ["reject", "real", "imag", "abs", "phase"])} = "reject"
    options.slicePoints (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.slicePoints, 2), mustBeFinite} = 64
end
palettes = struct( ...
    'bright', {{'#000000', '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'}}, ...
    'muted', {{'#000000', '#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD'}}, ...
    'vibrant', {{'#000000', '#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB'}}, ...
    'light', {{'#000000', '#EEDD88', '#99DDFF', '#44BB99', '#FFAABB', '#EECC66', '#77AADD', '#EE8866', '#DDDDDD'}}, ...
    'high_contrast', {{'#000000', '#004488', '#BB5566', '#DDAA33', '#FFFFFF'}}, ...
    'filbet', {{'#000000', '#CBCBFF', '#D04231', '#EF9F93', '#AE7AC5', '#5F94EE', '#9C640D', '#222F3D', '#8A8A8A'}});
if ~isfield(palettes, options.discreteColorScheme)
    error('mdemat:tools:visual:InvalidPalette', 'Unknown discrete color scheme.');
end
if isempty(options.lineScheme) || isempty(options.scatterScheme)
    error('mdemat:tools:visual:InvalidStyle', 'Line and scatter schemes must not be empty.');
end
style = struct('Palette', {palettes.(options.discreteColorScheme)}, ...
    'ColorMap', options.continuousColorScheme, ...
    'LineScheme', {options.lineScheme}, 'ScatterScheme', {options.scatterScheme}, ...
    'LineWidth', options.lineWidth, 'FontSize', options.fontSize, ...
    'AxisLabelFontSize', options.axisLabelFontSize, 'LegendFontSize', options.legendFontSize, ...
    'TitleFontSize', options.titleFontSize, 'LegendLocation', options.legendLocation, ...
    'ColorPadding', options.colorPadding, 'SymmetricColorLimits', options.symmetricColorLimits, ...
    'LogX', options.logX, 'LogY', options.logY, 'ComplexDisplay', options.complexDisplay, ...
    'SlicePoints', options.slicePoints);
sunset = [; ...
            54, 75, 154; ... % dark blue (strong negative)
            64, 90, 168; ... % slightly lighter blue
            74, 105, 183; ... % medium-dark blue
            84, 120, 198; ... % medium blue
            94, 135, 205; ... % light blue
            110, 150, 215; ... % lighter blue
            126, 166, 220; ... % very light blue
            142, 182, 225; ... % pale blue
            158, 198, 230; ... % very pale blue
            174, 214, 235; ... % blue-white
            190, 225, 238; ... % very pale blue
            206, 232, 240; ... % almost white blue
            222, 238, 242; ... % white-blue
            238, 244, 244; ... % near white
            245, 240, 220; ... % warm white (transition)
            252, 236, 196; ... % pale warm
            254, 228, 172; ... % very light orange
            254, 218, 148; ... % light orange
            254, 208, 124; ... % medium-light orange
            254, 198, 100; ... % medium orange
            253, 188, 88; ... % orange
            252, 178, 76; ... % orange-red
            251, 168, 64; ... % red-orange
            248, 150, 52; ... % red
            245, 132, 40; ... % darker red
            240, 114, 28; ... % dark red
            235, 96, 16; ... % darker red
            225, 78, 8; ... % very dark red
            215, 60, 0; ... % darkest red
            180, 40, 0; ... % maroon
            145, 20, 0; ... % dark maroon (strong positive)
            ] / 255;
burd = [; ...
            33, 102, 172; ... % dark blue (negative extreme)
            67, 147, 195; ... % blue
            146, 197, 222; ... % light blue
            209, 229, 240; ... % very light blue
            247, 247, 247; ... % white (neutral)
            253, 219, 199; ... % very light red
            244, 165, 130; ... % light red
            214, 96, 77; ... % red
            178, 24, 43; ... % dark red (positive extreme)
            ] / 255;
iridescent = [; ...
            254, 251, 233; ... % very light yellow (low values)
            252, 247, 213; ... % light yellow
            245, 243, 193; ... % yellow
            234, 240, 181; ... % yellow-green
            221, 236, 191; ... % light green
            208, 231, 202; ... % green
            194, 227, 210; ... % green-cyan
            181, 221, 216; ... % cyan
            168, 216, 220; ... % light cyan
            155, 210, 225; ... % cyan-blue
            141, 203, 228; ... % light blue
            129, 196, 231; ... % blue
            123, 188, 231; ... % medium blue
            126, 178, 228; ... % blue-purple
            136, 165, 221; ... % purple-blue
            147, 152, 210; ... % purple
            155, 138, 196; ... % light purple
            157, 125, 178; ... % purple-pink
            154, 112, 158; ... % pink-purple
            144, 99, 136; ... % dark purple
            128, 87, 112; ... % very dark purple
            104, 73, 87; ... % dark grey
            70, 53, 58; ... % very dark (high values)
            ] / 255;
maps = struct('sunset',sunset,'burd',burd,'iridescent',iridescent);
if isfield(maps, options.continuousColorScheme)
    style.ColorMap = maps.(options.continuousColorScheme);
end

end
