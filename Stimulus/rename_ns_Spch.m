

clc

src_path = 'D:\DataBases\Reverebration Stimuli\Project 1\ns_Spch_fc(4000Hz)\ns_Spch_fc(4000Hz)_Dist30_Rev2';
trg_path = 'D:\DataBases\Reverebration Stimuli\Project 1\ns_Spch_fc(4000Hz)\rep';

files = dir( src_path );

idx_files = arrayfun(@(S) ~S.isdir, files);
files = files(idx_files);

idx_wav = arrayfun(@(X) ~isempty(strfind(X.name,'.wav')), files, 'UniformOutput', 1);
files = files(idx_wav);

idx_brir = arrayfun(@(X) ~isempty(strfind(X.name,'ns_Spch')), files, 'UniformOutput', 1);
files = files(idx_brir);


%% Move the files
for kk = 1:length(files)
    fn = files(kk).name;
    
    % new filename
    expression = '4000Hz)';
    replace    = '4k)Hz';
    fn_new = regexprep(fn, expression, replace);
    if isempty(fn_new), continue; end
    %fn_new = fn;
    
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
    %params_new{1} = num2str( 1/1000*str2double(params{1}) );   % (cm)
    params_new{2} = num2str( 10*str2double(params{2}) );   % (cm)
    params_new{3} = num2str( 10*str2double(params{3}) );   % (%)
    fn_new = regexprep(fn_new, params, params_new);
    
    % Before saving
    fprintf('--> fn:     %s\n', fn);
    fprintf('--> fn_new: %s\n\n', fn_new);
        
    src = [files(kk).folder, filesep, fn];
    dst = [trg_path, filesep, fn_new];
    
    movefile(src, dst);
    
end



