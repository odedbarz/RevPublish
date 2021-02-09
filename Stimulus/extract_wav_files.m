function files = extract_wav_files(src_path)
%
%   function files = extract_wav_files(src_path)
%
%
% src_path: is a cell array of strings to folders that holds the WAV files
%           to extract.
%


if ischar(src_path), src_path = {src_path}; end

files = [];

for kk = 1:length(src_path)
    
    files_kk = dir( src_path{kk} );

    idx_files = arrayfun(@(S) ~S.isdir, files_kk);
    files_kk = files_kk(idx_files);
    if 0 == nnz(idx_files), continue; end

    idx_wav = arrayfun(@(X) ~isempty(strfind(X.name,'.wav')), files_kk, 'UniformOutput', 1);
    files_kk = files_kk(idx_wav);
    if 0 == nnz(idx_files), continue; end
    
    files_kk = arrayfun(@(S) [S.folder, filesep, S.name], files_kk, 'UniformOutput', 0);
    
    files = [files; files_kk];
    
end

