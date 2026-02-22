function renderVariable1d(axes, dataset, styleset, varName, offset)
subdata = prepareData(dataset, varName, 1);

if isempty(fieldnames(subdata))
    return;
end

[~] = render1d(axes, subdata, styleset, varName, offset);

end