function renderVariable2d(axes, dataset, styleset, varName, offset)
subdata = prepareData(dataset, varName, 2);

if isempty(fieldnames(subdata))
    return;
end

varAxes = render2d(axes, subdata, varName, offset);

renderSlices(axes, subdata, styleset, varName, offset, 2);

if length(varAxes) > 1
    syncColorLimits(varAxes, subdata);
end
end