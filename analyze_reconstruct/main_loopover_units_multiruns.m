%
% main_loopover_units_multiruns.m
%
% Description:
% This script loops over all loaded measurements and reconstruct the
% spectrogram using various size of database.
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
data_type   = 'SU';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(01-Nov-2021)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(08-Nov-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(10-Dec-2021)_bw(100)_fbands(30)_win(NaN)ms_spec(gammatone).mat'

fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data      = load(fn.load.fullfile);

spec_st   = data.spec_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));
n_units   = height(tbl_data);     % total available units
duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);


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
    aux.cprintf('UnterminatedStrings', '--> n_units     : %g (all available units)\n', n_units);
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



%%
if verbose
    fprintf('\n-> Starting training...\n');
    fprintf('========================\n');    
end

n_random_runs = 11;     % ************************* <<<<<<<<<<<< ========   
sort_type = 'SPK';
units = 10;             % ************************* <<<<<<<<<<<< ========   

sorted_list = find_best_unit_set(sort_type, 'fn_template', ...
    {fn.load.path, fn.load.file_template, data_type});
sort_type = 'SPK-RND'   % quick and dirty!

for q = 1    %1:n_drr 
    % The training (i.e., truth-level) DRR case
    train_drr = drr.sortby(q); %drr.dry;         
    %train_drr = [3, 4, 5];  '########## Training for more than one DRR #########'
    if verbose
        fprintf('--> TRAIN DRR: %d, %s\n', train_drr, drr.labels{train_drr});
    end    

    % Loop over UNITS
    for m = 1:n_random_runs
        fprintf('%d Starting a new random run...\n', m);

        %[sorted_list, ~] = find_best_unit_set(sort_type, 'N', n_units);
        % Randomize:
        rand_idx = randperm(length(sorted_list), units);
        H_sorted = data.H( :, 1:n_drr, sorted_list(rand_idx) );        
        
        % >> analyze_units;
        obj_list = cell(n_drr, n_splits);

        clear scores
        scores.CC  = nan(n_drr, n_splits);
        scores.mse = nan(n_drr, n_splits);
        scores.nmse= nan(n_drr, n_splits);
        scores2 = scores;
        
        % Get valid measurements
        y1 = squeeze( H_sorted(:, train_drr, :) );

        % Split the TRAINING set
        X1 = spec_st.Sft{train_drr};
        %X1 = [spec_st.Sft{train_drr}; spec_st.Sft_right{train_drr}];
        
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
                %test_drr  = k;
                test_drr = drr.ordered(k);                                               
                
                % Split for the TESTING data
                X2 = spec_st.Sft{test_drr};
                %X2 = [spec_st.Sft{test_drr}; spec_st.Sft_right{test_drr}];

                y2 = squeeze( H_sorted(:, test_drr, :) );                
                
                [~, X_test_kn, ~, y_test, ~] = train_test_split(X2, y2, ...
                    'n_splits', n_splits, ...
                    'split_time_idx', split_time_idx, ...
                    'test_grp', test_grp_number );

                % RECONSTRUCTION   
                % Predict the spectrogram
                obj.predict(y_test);
                

                % transform and keep as a structure for later analysis
                warning off
                obj_list{k,n}           = struct(obj);
                obj_list{k,n}           = rmfield(obj_list{k,n}, 'X_train');    % save memory
                obj_list{k,n}           = rmfield(obj_list{k,n}, 'r_train');    % save memory
                obj_list{k,n}.train_drr = train_drr;    % (1x1) index of the training DRR condition
                obj_list{k,n}.test_drr  = test_drr;     % (1x1) index of the testing DRR condition
                obj_list{k,n}.R_test    = [];           % convmtx of y_test

                % The reconstruction filter G is the same for all DRR
                % conditions, so save space on the HD
                if k > 1
                    obj_list{k,n}.G       = [];
                    obj_list{k,n}.G_sd    = [];
                    obj_list{k,n}.Crr     = [];
                    obj_list{k,n}.Crr_inv = [];
                    obj_list{k,n}.Crs     = [];
                end
                warning on

                % ### DEBUG ###
                % %{
                % Goodness-of-fit
                gof = goodness(X_test0, obj.X_est);
                scores.CC(k,n)  = gof.CC;
                scores.mse(k,n) = gof.mse;
                scores.nmse(k,n)= gof.nmse;
                
                % Goodness-of-fit
                gof2 = goodness(X_test_kn, obj.X_est);
                scores2.CC(k,n)  = gof2.CC;
                scores2.mse(k,n) = gof2.mse;
                scores2.nmse(k,n)= gof2.nmse;
                               
                %plot(1:5, scores.CC, 'o', 1:5, scores2.CC, 'x');
                %bar([mean(scores.CC,2), mean(scores2.CC,2)]);
                %}
            end
                
        end
    
    
    
    %% Save the reconstruction results
    % %{
        fprintf('SAVE the analysis data!\n');
        fn.save.path    = '../_data/Reconstruct/';
        fn.save.file    = sprintf('%d_reconstruct_%s_(%s)_units(%d)_bw(%g)ms_type(%s)_fbands(%d)_lags(%g)ms_cau(%d)',...
            m, data_type, date, units, binwidth, sort_type, n_bands, lags_ms, iscausal);
        fn.save.fullfile= fullfile( fn.save.path, fn.save.file );
                
        stim_st = data.stim_st;
        
        % Save the results for that 
        save(fn.save.fullfile, ... 
            '-v7.3', ...
            'splits', ...       saves the chunks\intervals\speakers
            'stim_st', ...      stimulus data
            'spec_st', ...      spectrogram's structue with all relevant data
            'obj_list', ...     cell array of reconstruction objetcs, saved as structures
            'tbl_data', ...     a table with all neurons in the data set            
            'fn', ...           filenames, including the data-set filename used here 
            'H_sorted',...      sorted units used for the analysis
            'sorted_list',...    the list of sorted unit used to create H_sorted from data.H
            'sort_type' ...
         );         
    %}
    end

end




