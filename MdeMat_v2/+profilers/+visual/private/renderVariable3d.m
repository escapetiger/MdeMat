function renderVariable3d(axes, dataset, styleset, varName, offset)
subdata = prepareData(dataset, varName, 3);

if isempty(fieldnames(subdata))
    return;
end
varAxes = render3d(axes, subdata, varName, offset);
renderSlices(axes, subdata, styleset, varName, offset, 3);
if length(varAxes) > 1
    syncColorLimits(varAxes, subdata);
end
end