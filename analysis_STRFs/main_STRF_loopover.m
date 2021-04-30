%
% main_STRF_loopover.m
%
% Description:
% This script loops over all loaded measurements and reconstruct the
% spectrogram using various size of database.
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
data            = load(fn.load.fullfile);
n_units         = size(data.H,3); 
spec_st         = data.spec_st;
tbl_data        = data.(sprintf('tbl_%s', data_type));
duration_sec    = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_stats.m] Loading file:\n\t...<%s>\n', fn.load.file);





%%
iscausal    = 1;                    % use reconstructed causal filters? 
lags_ms     = 30;                   % (ms) maximum lags
binwidth    = spec_st.binwidth;     % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 
jk_flag     = 1;     
algo_type   = 'regression';    % {'regression', 'asd'}


if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> n_units     : %d\n', n_units);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    
    aux.cprintf('UnterminatedStrings', '    STRFs:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife: %d\n', jk_flag);
end




%%
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));



 
%% Create & initialize the STRF object
obj = strf_c(binwidth,...
    'f', spec_st.f,...
    'iscausal', iscausal, ...
    'lags_ms', lags_ms,...
    'jk_flag', jk_flag); 



%% Start looping over all cases
for q = 1 %1:drr.n_drr 
    % The training (i.e., truth-level) DRR case
    train_drr = drr.sortby(q); %drr.dry;      
    if verbose
        fprintf('--> TRAIN DRR: %d, %s\n', train_drr, drr.labels{train_drr});
    end

    % Select M_UNITS to reconstruct
    H_train = squeeze(data.H(:, train_drr, :));


    % Save all STRFs stuff in one table 
    tbl_strf = array2table(cell(n_units,4), ...
        'VariableNames', {'strf', 'strf_std', 'strf_speaker', 'y_strf_dry'} );    
    
    % Add best frequency
    bf = zeros(n_units,1);
    bf_std = zeros(n_units,1);
    tbl_strf = [tbl_strf, table(bf, bf_std)];
    clear bf bf_std
    
    % Save the CC
    CCest_dry = cell(n_units,1);
    tbl_strf = [tbl_strf, array2table(CCest_dry)];        
    clear CCest_dry
    
    CCest_drr = cell(n_units,1);
    tbl_strf = [tbl_strf, array2table(CCest_drr)];
    clear CCest_drr
    
    % Add row number to the beginning of the table
    index = (1:n_units)';
    tbl_strf = [array2table(index), tbl_strf];
    clear index
    
    
    %% Loop over UNITS
    for m = 1:n_units
        if verbose && 0 == rem(m,2)
            fprintf('\n--> unit #%d (out of %d)\n', m, n_units);
        end

        BF_speaker = nan(1, n_splits);   % Best frequency for each split
        strf_speaker = nan(n_bands, obj.n_lags, n_splits);
        
        % Get valid measurements
        y1 = data.H(:, train_drr, m);
        
        % Split the TRAINING set
        X1 = spec_st.Sft{train_drr};

        CCnk = nan(n_splits, n_drr);        % est-to-dry
        CC2nk = nan(n_splits, n_drr);       % est-to-drr

        y_strf_dry = zeros(length(y1), drr.n_drr);
        
        
        %% Loop overSPLITS
        for n = 1:n_splits
            test_grp_number = n;
            if verbose && 0 == rem(m,2) && 0 == rem(n,4)
                fprintf('----> speaker #: (%d/%d)\n', test_grp_number, n_splits);
            end

            [X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
                'split_time_idx', split_time_idx,...
                'test_grp', test_grp_number ...
            );
            
            % Fit the model
            [strf_speaker(:,:,n), ~] = obj.fit(X_train, y_train,...
                'fignum', [],...
                'verbose', []);
            
            % Calculate the best frequency
            [~, idx_bf_n] = max( max(squeeze(strf_speaker(:,:,n)),[],2) );     % (Hz)
            BF_speaker(n) = spec_st.f(idx_bf_n);

            % Loop over DRRs
            for k = 1:drr.n_drr    
                % Choose the DRR case for the testing signals
                test_drr = drr.ordered(k);

                % Split for the TESTING data
                X2 = spec_st.Sft{test_drr};
                y2 = data.H(:, test_drr, m);
                
                [~, X_test, ~, y_test, ~] = train_test_split(X2, y2, ...
                    ...'n_splits', n_splits, ...   
                    'split_time_idx', split_time_idx,...
                    'test_grp', test_grp_number );

                % *** STRF ***
                % Predict the spectrogram
                r_est = obj.predict(X_test);
                gof  = goodness(y_test0, r_est);     % est-to-dry
                gof2 = goodness(y_test, r_est);      % est-to-drr
                
                % Goodness-of-fit
                CCnk(n,k)  = gof.CC;       % est-to-dry
                CC2nk(n,k) = gof2.CC;      % est-to-drr
                
                % Aggregate all reference (aka dry) responses
                idx_split = split_time_idx(test_grp_number,1):split_time_idx(test_grp_number,2);
                y_strf_dry(idx_split, k) = r_est;
            end

        end
    
        tbl_strf.CCest_dry{m}   = CCnk;         % est-to-dry     
        tbl_strf.CCest_drr{m}   = CC2nk;        % est-to-drr   
        tbl_strf.bf(m)          = mean(BF_speaker);          % (Hz)
        tbl_strf.bf_std(m)      = std(BF_speaker);           % (Hz)
        tbl_strf.strf{m}        = mean(strf_speaker, 3);
        tbl_strf.strf_speaker{m}= strf_speaker;
        tbl_strf.strf_std{m}    = std(strf_speaker, [], 3);
        tbl_strf.y_strf_dry{m}  = y_strf_dry;
    end
    
    
    
    %% Save the reconstruction results
    % %{
        obj = strf_c(binwidth,...
            'f', spec_st.f,...
            'iscausal', iscausal, ...
            'lags_ms', lags_ms,...
            'jk_flag', jk_flag ); 
        strf_st = struct(obj);
    
        fprintf('SAVE the STRF''s analyzed data!\n');
        fn.save.path    = '../_data/Analysis/';
        fn.save.file    = sprintf('STRF_%s_(%s)_units(%d)_bw(%g)ms_algo(%s)_fbands(%d)_splits(%d)_lags(%g)ms_cau(%d)_trainDRR(%s)',...
            data_type, date, n_units, binwidth, algo_type, n_bands, n_splits, lags_ms, iscausal, num2str(train_drr, '%d '));
        fn.save.fullfile= fullfile( fn.save.path, fn.save.file );
                
        % Save the results for that 
        save(fn.save.fullfile, '-v7.3', ...
            'splits', ...       saves the chunks\intervals\speakers
            'spec_st', ...      spectrogram's structue with all relevant data
            'fn', ...           filenames, including the data-set filename used here 
            'strf_st',...       STRF_C object
            'tbl_strf'  ...     (table) all STRFs and other data (e.g., BFs)
            ); 
        
    %}
    

end




