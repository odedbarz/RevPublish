%
% main_loopover_stats.m
%
% Description:
% Loop over (randomly) selected units and saves the result.
%
% Notes:
% * Before running this file you need to run main_aggregate_MUA_data.m.
% This will create a data file at ./_data of the measurements to
% reconstuct.
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');
   
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'MUA';       % {'SU', MUA'}
fn.load.path= '../_data';
data_type   = upper(data_type);

switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn.load.file = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        n_units = 50; %[10, 25, 50, 103];

    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        n_units = 50; %[10, 25, 50, 103, 241];

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data = load(fn.load.fullfile);

spec_st         = data.spec_st;
tbl_data        = data.(sprintf('tbl_%s', data_type));
len_unit_list   = length(n_units);
duration_sec    = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_stats.m] Loading file:\n\t...<%s>\n', fn.load.file);



%% Reconstruction parameters
iscausal       = 0;                    % use reconstructed causal filters?
lags_ms        = 30;                   % (ms) maximum lags
binwidth       = spec_st.binwidth;     % (ms)
n_bands        = spec_st.n_bands;
win_size_ms    = spec_st.win_size_ms; 
jackknife_flag = 1;     
algo_type      = 'regression';    % {'regression', 'asd', 'svd'}

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type       : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> n_units         : [%s]\n', num2str(n_units));
    aux.cprintf('UnterminatedStrings', '--> height(tbl_data): %g\n', height(tbl_data));
    aux.cprintf('UnterminatedStrings', '--> duration_sec    : %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality       : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms         : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth        : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands         : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms     : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife    : %d\n', jackknife_flag);
end



%%
%{
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));


'### DEBUG ###'
splits_range = 1:6 %1:12 %1:6 %7:12;
n_splits = length(splits_range);
split_time_idx = split_time_idx(splits_range,:);
split_time_idx = split_time_idx - split_time_idx(1,1) + 1;
'### DEBUG ###'
%}


%%
if verbose
    fprintf('\n-> Starting training...\n');
    fprintf('========================\n');    
end

% Number of repetitions for later statisic analysis
n_rep = 30
splits_range = 0 + (1:11),  '######'
n_splits = length(splits_range);

aux.cprintf('UnterminatedStrings', '--> n_rep   : %d\n', n_rep);

Hdry = squeeze( data.H(:,drr.ordered(1),:) ); 

CCstats = nan(n_rep, n_drr, n_splits);
% CCunits = nan(n_units, n_rep);

'### DEBUG ###'
clear CCunits split_time_idx_reps
dummy = load(fullfile(load.path_to_data('Stats'), 'CCunits_(09-Jun-2021)_for_valid_stats_test.mat'));   % --> CCunits
CCunits = dummy.CCunits;
split_time_idx_reps = dummy.split_time_idx_reps;




%% !!!IMPORTANT!!! Use ONE time to create a reference set for the analysis
%{
n_all_splits = 12;

% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_all_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));

CCunits = nan(n_units, n_rep);
splits_range = nan(n_all_splits, n_rep);
split_time_idx_reps = nan(n_all_splits, 2, n_rep);
split_time_diff = diff(split_time_idx,[],2);
for q = 1:n_rep 
    CCunits(:,q) = randperm(size(data.H,3), n_units);  
    
    splits_range(:,q) = randperm(n_all_splits);
    split_time_idx_q = split_time_diff(splits_range(:,q));
    split_time_idx_q = [[0; 1; split_time_idx_q(1:end-1)+2], [0; split_time_idx_q+1]];
    split_time_idx_q = cumsum(split_time_idx_q);
    split_time_idx_q = split_time_idx_q(2:end,:);
    
    split_time_idx_reps(:,:,q) = split_time_idx_q;
end

% !!! SAVE !!!
path_ccunits = load.path_to_data('Stats');
save(fullfile(path_ccunits, 'CCunits_(00-XXX-2021)_for_valid_stats_test.mat'),...
    'CCunits', ...
    'splits_range', 'split_time_idx_reps');
%}


%%
for q = 1:n_rep    
    % Select a new set of units for each iterations
    %sorted_list = find_best_unit_set('RND', 'Y', Hdry);
    %rng('default');     % Generate random numbers that are repeatable
    %sorted_list = randperm(size(data.H,3), n_units);    
    %CCunits(:,q) = sorted_list(1:n_units);

    split_time_idx = split_time_idx_reps(splits_range,:,q);
    split_time_idx = split_time_idx - split_time_idx(1,1) + 1;
    
    % Select M_UNITS to reconstruct        
    H_sorted = data.H( :, 1:n_drr, CCunits(:,q) );      
    
    % The training (i.e., truth-level) DRR case
    train_drr = drr.dry;         
    %train_drr = drr.ordered(end);

    if verbose
        fprintf('--> TRAIN DRR: %d, %s\n', train_drr, drr.labels{train_drr});
    end

    % >> analyze_units;
    obj_list = cell(n_drr, n_splits);


    % Loop over SPLITS
    for n = 1:n_splits
        % choose the testing chunk\speaker out of the stimulus
        test_grp_number = n;    
        if verbose
            fprintf('----> splits #: (%d/%d)\n', test_grp_number, n_splits);
        end

        % Enable to train on more than 1 DRR            
        % %{
        assert(1 == length(train_drr), '--> Use this Option for only ONE train_drr!');

        % Get valid measurements
        y1 = squeeze( H_sorted(:, train_drr, :) );

        % Split the TRAINING set
        X1 = spec_st.Sft{train_drr};
        
    '### DEBUG ###'
    y1 = y1(split_time_idx(1,1):split_time_idx(end,2),:);
    X1 = X1(:, split_time_idx(1,1):split_time_idx(end,2));
    '### DEBUG ###'

        [X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
           'n_splits', n_splits, ...
           'split_time_idx', split_time_idx, ...
           'test_grp', test_grp_number ...
        );
        %}


        % Extract the reconstruction filters
        % Initialize the reconstruction object
        obj = reconstruct_c(binwidth,...
            'f', spec_st.f, ...
            'iscausal', iscausal, ...   
            'algo_type', algo_type, ... {'regression', 'asd'}                
            'lags_ms', lags_ms ); 

        % Fit the model
        obj.fit(X_train, y_train, ...
            'jk_flag', jackknife_flag, ...
            'n_splits', n_splits, ...
            'fignum', []);


        % Loop over DRRs
        for k = 1:n_drr    
            % Choose the DRR case for the testing signals
            test_drr  = drr.ordered(k);
            if 0 % verbose
                fprintf('-----> TEST DRR: %s\t(%d/%d)\n',  drr.labels{test_drr}, test_drr, n_drr);
            end

            % Split for the TESTING data
            X2 = spec_st.Sft{test_drr};
            y2 = squeeze( H_sorted(:, test_drr, :) );
            
    '### DEBUG ###'
    y2 = y2(split_time_idx(1,1):split_time_idx(end,2),:);
    X2 = X2(:, split_time_idx(1,1):split_time_idx(end,2));
    '### DEBUG ###'
            
            [~, X_test_kn, ~, y_test, ~] = train_test_split(X2, y2, ...
                'n_splits', n_splits, ...
                'split_time_idx', split_time_idx, ...
                'test_grp', test_grp_number );

            % RECONSTRUCTION -- Predict the spectrogram
            obj.predict(y_test);

            % Goodness-of-fit
            gof = goodness(X_test_kn, obj.X_est);
            CCstats(q,k,n) = gof.CC;
            
        end


    end
    

end

% DRRs in CCstats are ** ORDERED **
drr_labels = drr.labels(drr.ordered);


%% Save the reconstruction results
% %{
    fprintf('SAVE the analysis data!\n');
    fn.save.path    = '../_data/Stats/';
    fn.save.file    = sprintf('CCstats_TEST(all)_%s_(%s)_units(%d)_bw(%g)ms_algo(%s)_fbands(%d)_splits(%d)_lags(%g)ms_cau(%d)_trainDRR(%s)',...
        data_type, date, n_units, binwidth, algo_type, n_bands, n_splits, lags_ms, iscausal, num2str(train_drr, '%d '));
    fn.save.fullfile= fullfile( fn.save.path, fn.save.file );

    stim_st = data.stim_st;

    % Save the results for that 
    save(fn.save.fullfile, ... 
        'drr_labels', ...   info: DRRs in CCstats are ** ORDERED **
        'CCstats', ...      saves the chunks\intervals\speakers
        'CCunits',...       the list of sorted unit used to create H_sorted from data.H
        'stim_st', ...      stimulus data
        'spec_st', ...      spectrogram's structue with all relevant data
        'tbl_data', ...     a table with all neurons in the data set            
        'fn', ...           filenames, including the data-set filename used here 
        'H_sorted'...       sorted units used for the analysis
        ); 

%}



