% create_dictionary.m
%
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');




%% Load spectrogram data
% %{
data_type = 'MUA';       % {'SU', MUA'}
path_stim = '../_data';
data_type = upper(data_type);
switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        filename = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        filename = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn_2_load = fullfile( path_stim, filename );
dummy = load(fn_2_load, 'spec_st');
spec_st = dummy.spec_st;

n_bands = spec_st.n_bands;


%% Load signals to build the dictionary
% Path to the Impale's data
%path_root_mat = load.path_to_data('Impale_data');

path_signals = fullfile( load.path_to_data('data'), 'Dictionary', 'signals' );
files = dir(path_signals);

valid_file = false(1, length(files));
for k = 1:length(files)
    if files(k).isdir
        continue;
    end
    
    [~, ~, ext] = fileparts( files(k).name );
    if ~strcmpi(ext, '.wav')
        continue;
    end
    
    valid_file(k) = true;
        
end
files = files( valid_file );





%%
Strain = [];

for k = 1:length(files)
    fn = files(k).name;
    [ywav, fs] = audioread(fullfile(path_signals, fn));
    wav_info   = audioinfo(fullfile(path_signals, fn));

    [Sk, ~] = spec.spectrogram(ywav, fs, ...
        'n_bands', spec_st.n_bands,...
        'lowfreq', spec_st.lowfreq,...
        'highfreq', spec_st.highfreq,...
        'overlap_ratio', spec_st.overlap_ratio,...
        'binwidth', spec_st.binwidth,...
        'win_size_ms', spec_st.win_size_ms, ...
        'nw', spec_st.nw,...                only for spectrogram_type== MULTITAPER
        'f_scale', spec_st.f_scale,...
        'db_floor', spec_st.db_floor, ...  % (dB)
        'method', spec_st.method, ...
        'fignum', [] ...
     );
    %title(aux.mName2latex(fn));
    %pause(0.001);
    
    % Concatenate all the signals
    Strain = [Strain, Sk];    
    
end




%% Create the dictionary
lags_ms     = 30;                       % (ms) the time duration of each patch 
binwidth    = spec_st.binwidth;         % (ms) 
n_lags      = fix(lags_ms/binwidth);    % (smp) number of samples of each patch (along the time duration)
n_train     = length(Strain);           % (smp) 
n_features  = 13*13;


dict_type = 'NNMF';                 % {'FASTICA', 'RICA', 'PCA', 'DL', 'NNMF'}
fprintf('--> Calculating dictionary...\n');

[D, dict_st] = compute_bases(Strain, n_lags,...
    'dict_type', dict_type, ...
    'n_features', n_features, ...
    'verbose', 0, ...
    'fignum', 1);
    
fprintf('--> Finished\n');

fprintf('\nDictionary:\n');
fprintf('-> lags_ms     : %g ms\n', lags_ms);
fprintf('-> binwidth    : %g ms\n', binwidth);
fprintf('-> n_lags      : %d\n', n_lags);
fprintf('-> n_train     : %d\n', n_train);
fprintf('-> dict_type   : %s\n', dict_type);




%% Recorded response to DRY spectrograms:
drr     = get_DRR_list_and_indices;
n_drr   = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);

idx_test= idx_dry;
Stest = spec_st.Sft{idx_test};
Stest_col = im2col(Stest, [n_bands, n_lags], 'distinct');

Xtest = Stest_col(:, randperm(n_features));
Xtest = Xtest(:, sum(abs(Xtest),1)>1e-6 );  	% avoid zero patches
n_samples = size(Xtest, 2);

% recorded response
RS = D' * Xtest/size(D,1);






%% Save the dictionary
fprintf('saving the dictionary...\n');

filename = sprintf('DIC_%s_(%s)_type(%s)_lags(%g)ms_bw(%g)_fbands(%d)_win(%g)ms_spec(%s)', ...
    data_type, date, dict_type, lags_ms, binwidth, spec_st.n_bands, spec_st.method );
fn_2_save = fullfile( load.path_to_data('data'), 'Dictionary', filename);

fprintf('\nSaving data at:\n');
disp(['--> ', fn])
save(fn_2_save, 'D', 'dict_st', 'RS', 'spec_st');

fprintf('Done!\n')





