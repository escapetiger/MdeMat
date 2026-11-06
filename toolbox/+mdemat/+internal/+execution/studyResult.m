function report = studyResult(scheme, state)
%STUDYRESULT Return the accepted output of a retained research run.
report = struct('schemaVersion', 1, 'scheme', class(scheme), ...
    'completed', true, 'diagnostics', state.History);
if isprop(scheme, 'TDisc')
    report.finalTime = scheme.TDisc.Timeline.Now;
    report.completed = scheme.TDisc.Timeline.IsExhausted;
end
end
