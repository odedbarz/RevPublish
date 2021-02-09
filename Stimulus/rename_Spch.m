

clc

src_path = 'D:\DataBases\Reverebration Stimuli\Project 1\Spch\BRIR_Dist30_Rev8';
trg_path = 'D:\DataBases\Reverebration Stimuli\Project 1\Spch\reps';

files = dir( src_path );

idx_files = arrayfun(@(S) ~S.isdir, files);
files = files(idx_files);

idx_wav = arrayfun(@(X) ~isempty(strfind(X.name,'.wav')), files, 'UniformOutput', 1);
files = files(idx_wav);

idx_brir = arrayfun(@(X) ~isempty(strfind(X.name,'BRIR')), files, 'UniformOutput', 1);
files = files(idx_brir);


%% Move the files
for kk = 1:length(files)
    fn = files(kk).name;
    
    % new filename
    expression = 'BRIR';
    replace    = 'Spch';
    fn_new = regexprep(fn, expression, replace);
    if isempty(fn_new), continue; end
    
    expression = '-Dist';
    replace    = '-dst';
    fn_new = regexprep(fn_new, expression, replace);
    if isempty(fn_new), error('!'); end

    expression = '-Rev';
    replace    = '-rev';
    fn_new = regexprep(fn_new, expression, replace);
    if isempty(fn_new), error('!'); end

    expression = '-Rep';
    replace    = '-rep';
    fn_new = regexprep(fn_new, expression, replace);
    if isempty(fn_new), error('!'); end
    
    % Change the units
    params = regexp(fn_new,'\d*','Match');
    params_new = params;
    params_new{1} = num2str( 10*str2double(params{1}) );   % (cm)
    params_new{2} = num2str( 10*str2double(params{2}) );   % (%)
    fn_new = regexprep(fn_new, params, params_new);
    
    % Before saving
    fprintf('--> fn:     %s\n', fn);
    fprintf('--> fn_new: %s\n\n', fn_new);
        
    src = [files(kk).folder, filesep, fn];
    dst = [trg_path, filesep, fn_new];
    
    movefile(src, dst);
    
end



