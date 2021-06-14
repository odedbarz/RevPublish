%
% main_loopover_units.m
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
switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        %fn.load.file = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        unit_list = [10, 25, 50, 103];

    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        %fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';  
        unit_list = [10, 25, 50, 103]; %[10, 25, 50, 103, 241];

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data      = load(fn.load.fullfile);

spec_st   = data.spec_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));
n_units   = height(tbl_data);     % total available units
len_unit_list = length(unit_list);
duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);



%% Order the units
% [sorted_list, tbl_BF] = find_best_unit_set('CC', 'fn', fn.load.fullfile);

sort_type = 'CC';  % {'CC', 'RND', 'SVD', 'FILE', 'SPK', 'NOSPK'}
% Hdry = squeeze( data.H(:,drr.ordered(1),:) ); 
% sorted_list = find_best_unit_set(sort_type, 'Y', Hdry);                 % {'RND'}
% sorted_list = find_best_unit_set(sort_type, 'Y', Hdry, 'n_svd', 10);    % {'RND'}
[sorted_list, tbl_BFcc] = find_best_unit_set(sort_type, 'fn', fn.load.fullfile);  % {'CC'}
% sorted_list = find_best_unit_set(sort_type, 'fn', 'unit_list_MUA_drr(5).mat');  % {'FILE'}
% find_best_unit_set('FILE', 'fn', 'idx_MUA_good_sorted_unit_thr(0-7).mat'); % {'FILE'}
% sorted_list = find_best_unit_set(sort_type, 'fn_template', ...
%     {fn.load.path, fn.load.file_template, data_type});    % {'SPK', 'NOSPK'}


%% Get best measurement (CC) for each frequency group
% %{
    T = tbl_BFcc;
%         dummy = load('/Users/ob993/Dropbox (Partners HealthCare)/codeOnCloud/RevPublish/_data/Analysis/STRF_SU_(22-Apr-2021)_units(103)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat');
%         tbl_strf = dummy.tbl_strf;
%         tbl_strf.BF = tbl_strf.bf;
%         tbl_strf.neuron = tbl_strf.index;
%         tbl_strf.R = log10(1./(tbl_strf.bf_std + eps));
%         T = tbl_strf;
    n_grps = 10;
    frange = linspace(min(spec_st.f), max(spec_st.f), n_grps);

    % divide units into frequency groups
    bf_grps = cell(1, n_grps);
    I1 = T.BF >= frange;
    %I2 = T.BF <= frange;
    %Igrp = I1 .* circshift(I2,-1,2);
    Igrp = sum(I1, 2);
    for k = 1:n_grps
        idx_neurons = find(k==Igrp);
        
        % Order such that best CC are at the top of the list
        cc_neurons_k = T.R(idx_neurons);
        [~, idx_neurons_sort_k] = sort(cc_neurons_k, 'descend'); 
        bf_grps{k} = idx_neurons(idx_neurons_sort_k);
        
        %'### DEBUG : DONT SORT ###'
        %bf_grps{k} = idx_neurons
    end
    assert(sum( arrayfun(@(X) length(bf_grps{X}), 1:n_grps) ) == n_units);


    % Sort the neurons; assuming the all groups are non-empty at the beginning
    sorted_list = nan(1, n_units);
    klist = 1;     % the sorted_list counter
    J = 1;  % counter
    while klist <= n_units
        idx_grp = mod(J-1, n_grps)+1;        
        neuron = bf_grps{idx_grp};
        if ~isempty(neuron)
            sorted_list(klist) = neuron(1);
            klist = klist + 1;
            bf_grps{idx_grp} = neuron(2:end);            
        end
        J = J + 1;        
    end
    assert(length(sorted_list) == n_units);
    
%}


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


% % Use these 2 lines to override the split by speakers:
% split_time_idx = [];
% n_splits = 200




%%
if verbose
    fprintf('\n-> Starting training...\n');
    fprintf('========================\n');    
end

    

for q = 1 %1:n_drr 
    % The training (i.e., truth-level) DRR case
    train_drr = drr.sortby(q); %drr.dry;         
    %train_drr = [3, 4, 5];  '########## Training for more than one DRR #########'
    if verbose
        fprintf('--> TRAIN DRR: %d, %s\n', train_drr, drr.labels{train_drr});
    end    
    
    % Loop over UNITS
    for m = 1:len_unit_list
        % Set the # of neurons for the reconstruction
        m_units = unit_list(m);
        if verbose
        end

        % Select M_UNITS to reconstruct        
        H_sorted = data.H( :, 1:n_drr, sorted_list(1:m_units) );        
        
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
        fn.save.file    = sprintf('reconstruct_%s_(%s)_units(%d)_bw(%g)ms_algo(%s)_fbands(%d)_splits(%d)_lags(%g)ms_cau(%d)_trainDRR(%s)',...
            data_type, date, m_units, binwidth, algo_type, n_bands, n_splits, lags_ms, iscausal, num2str(train_drr, '%d '));
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




