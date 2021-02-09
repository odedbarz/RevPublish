function Twords = extract_TIMIT_meta( wav_file_list, path2TIMIT )
%
%   function Twords = extract_TIMIT_meta( wav_file_list, path2TIMIT )
%
% Description:
% Extract information (sentences, phonetics , etc.) from the TIMIT corpus
% meta-data fiels.
%

% Seperate between directories & filenames
[path2file, fn, file_ext] = cellfun(@(S) fileparts(S), wav_file_list, 'UniformOutput', 0);

% filename must include extensions
fn = cellfun(@(S1,S2) [S1, S2], fn, file_ext, 'UniformOutput', 0);  

Twords = table(fn, 'VariableName', {'fn'});

% PHN
Twords.phn = cell(size(Twords,1),1);            % *.PHN files
Twords.phn_start = cell(size(Twords,1),1);      % *.PHN files
Twords.phn_end = cell(size(Twords,1),1);        % *.PHN files
%WRD
Twords.wrd = cell(size(Twords,1),1);            % *.WRD files
Twords.wrd_start = cell(size(Twords,1),1);      % *.WRD files
Twords.wrd_end = cell(size(Twords,1),1);        % *.WRD files
% TXT
Twords.txt = cell(size(Twords,1),1);            % *.TXT files

% Save the path to the files
Twords.path = path2file;


% Parse the information in the TIMIT database
for kk = 1:length(wav_file_list)
    % Get the path to read the files from
    [~, fn_kk, ~] = fileparts( wav_file_list{kk} );    % remove the path to the filename (if exist)
    %fn_kk = [fn_kk, ext];
    path_in_TIMIT_subfolders = regexp(fn_kk, '_', 'split');    % get the TIMIT sub-folders
    
    path2files = cellfun(@(S) [S, filesep], path_in_TIMIT_subfolders(1:end-1), 'UniformOutput', 0);
    path2files = [path2TIMIT, filesep, [path2files{:}]];
    
    file2load = path_in_TIMIT_subfolders{end};
    
    % Read the TXT file
    fn_txt = [path2files, file2load, '.txt'];
    assert(~isempty(dir(fn_txt)), '--> ERROR: couldn''t find this file!');
    txt_k = textread(fn_txt, '%s', 'whitespace', ''); 
    Twords.txt{kk} = txt_k{1};
    
    % Read the PHN file
    fn_phn = [path2files, filesep, file2load, '.phn'];
    assert(~isempty(dir(fn_phn)), '--> ERROR: couldn''t find this file!');
    [Twords.phn_start{kk}, Twords.phn_end{kk}, Twords.phn{kk}] = textread(fn_phn, '%d %d %s'); 
    
    % Read the WRD file
    fn_wrd = [path2files, filesep, file2load, '.wrd'];
    assert(~isempty(dir(fn_wrd)), '--> ERROR: couldn''t find this file!');
    [Twords.wrd_start{kk}, Twords.wrd_end{kk}, Twords.wrd{kk}] = textread(fn_wrd, '%d %d %s'); 
    
end
