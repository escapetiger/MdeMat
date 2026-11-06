function sample = readStudySamples(filename)
%READSTUDYSAMPLES Read physical points and sampled values as primitive data.
sample = mdemat.tools.io.readRecord(filename);
fields = {'Points'; 'Values'; 'schemaVersion'; 'type'};
if ~isstruct(sample) || ~isscalar(sample) || ~isequal(sort(fieldnames(sample)), sort(fields)) ...
        || ~isequal(sample.type, "mdemat.sampledStudy") || ~isequal(sample.schemaVersion, 1) ...
        || ~isnumeric(sample.Points) || ~ismatrix(sample.Points) ...
        || ~isstruct(sample.Values) || ~isscalar(sample.Values)
    error('mdemat:study:InvalidSamples', 'Expected a versioned sampled study record.');
end
names = fieldnames(sample.Values);
for index = 1:numel(names)
    value = sample.Values.(names{index});
    if ~isnumeric(value) || ~ismatrix(value) || size(value, 1) ~= size(sample.Points, 2)
        error('mdemat:study:InvalidSamples', 'Each field requires one row per physical point.');
    end
end
end
