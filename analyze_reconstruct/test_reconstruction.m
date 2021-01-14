%
% test_reconstruction.m
%
% Description:
% Run reconstruction on ONE selected measurement
%
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
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data        = load(fn.load.fullfile);

spec_st   = data.spec_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));
n_units   = height(tbl_data);     % total available units
% unit_list = [10, 25, 50, 103, 150, n_units];
unit_list = 50,  '########### DEBUG ##############'


% Don't overflaw number of units in the database
unit_list = unique(min(height(tbl_data), unit_list));
len_unit_list = length(unit_list);


duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');


aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);



%% Order the units
% Get a list of VALID units to reconstruct from
% valid_units = squeeze( data.H(:,train_drr(1),:) );  % for VALID units
[sorted_list, tbl_BF] = find_best_unit_set('CC', 'fn', fn.load.fullfile);



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
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> unit_list   : [%s]\n', num2str(unit_list));
    aux.cprintf('UnterminatedStrings', '--> n_units     : %g\n', n_units);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife: %d\n', jackknife_flag);
end



%%
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));


% split_time_idx = [];
% n_splits = 200




%%
fprintf('\n-> Starting training...\n');


% The training (i.e., truth-level) DRR case
train_drr = drr.ordered(1); %drr.dry;         
fprintf('--> TRAIN DRR: %d, %s\n', train_drr, drr.labels{train_drr});
    
    
% Set the # of neurons for the reconstruction
m_units = unit_list(1);

% Select M_UNITS to reconstruct        
H_units = data.H( :, 1:n_drr, sorted_list(1:m_units) );        

        

%% >> analyze_units;
obj_list = cell(n_drr, n_splits);

        

%% Reconstruction
for n = 4 %1:n_splits
    % choose the testing chunk\speaker out of the stimulus
    test_grp_number = n;    
    fprintf('----> splits #: (%d/%d)\n', test_grp_number, n_splits);

    % Enable to train on more than 1 DRR            
    assert(1 == length(train_drr), '--> Use this Option for only ONE train_drr!');

    % Get valid measurements
    y1 = squeeze( H_units(:, train_drr, :) );

    % Split the TRAINING set
    X1 = spec_st.Sft{train_drr};

    [X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
       'n_splits', n_splits, ...
       'split_time_idx', split_time_idx, ...
       'test_grp', test_grp_number ...
    );



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
    for k = 4  %1:n_drr    
        % Choose the DRR case for the testing signals
        test_drr  = k;
        if 0 % verbose
            fprintf('-----> TEST DRR: %s\t(%d/%d)\n',  drr.labels{test_drr}, test_drr, n_drr);
        end

        % Split for the TESTING data
        X2 = spec_st.Sft{test_drr};
        y2 = squeeze( H_units(:, test_drr, :) );
        [~, X_test_kn, ~, y_test, ~] = train_test_split(X2, y2, ...
            'n_splits', n_splits, ...
            'split_time_idx', split_time_idx, ...
            'test_grp', test_grp_number );

        % RECONSTRUCTION   
        % Predict the spectrogram
        obj.predict(y_test);
        
        X_est = obj.X_est;  
        
        gof = goodness(X_test_kn, obj.X_est);
    end


end


%%
fprintf('\nResults -- Goodness-of-fit:\n');
disp(gof)


