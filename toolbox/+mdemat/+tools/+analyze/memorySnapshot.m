function report = memorySnapshot()
%MEMORYSNAPSHOT Report native Windows MATLAB process memory in bytes.
%   Bytes is memory().MemUsedMATLAB from this MATLAB process. It is one
%   snapshot after this function starts, not a high-water mark or allocation
%   count. PeakBytes remains NaN and PeakMeasured is false. On platforms where
%   MATLAB's documented memory API is unavailable, Available is false and
%   Bytes is NaN with a Reason. WSL launcher memory is never sampled.
report = struct('Available', false, 'Bytes', NaN, 'Units', 'bytes', ...
    'Kind', 'snapshot', 'PeakMeasured', false, 'PeakBytes', NaN, ...
    'Process', 'current MATLAB', 'Computer', computer, 'Reason', '');
if ~ispc
    report.Reason = 'The documented MATLAB memory API is supported on Windows.';
    return;
end
try
    usage = memory;
    report.Bytes = double(usage.MemUsedMATLAB);
    report.Available = true;
catch exception
    report.Reason = sprintf('%s: %s', exception.identifier, exception.message);
end
end
