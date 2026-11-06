function writeBuildEvidence(file, report)
%WRITEBUILDEVIDENCE Write an explicit build report as UTF-8 JSON.
directory = fileparts(file);
if ~isfolder(directory), mkdir(directory); end
stream = fopen(file, 'w', 'n', 'UTF-8');
if stream < 0, error('mdemat:build:EvidenceWrite', 'Cannot write %s.', file); end
guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
fprintf(stream, '%s\n', jsonencode(report, PrettyPrint=true));
end
