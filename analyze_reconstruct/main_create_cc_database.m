%
% main_create_cc_database.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions
% and saves the database into ANALYZED_CC_XX mat file.
% 
%
%

close all;
clc;
clear all;


fignum = 10;
verbose = 1;

setup_environment('../');



%% Load data
data_type   = 'SU'       % {'SU', MUA'}
fn_path     = load.path_to_data('Reconstruct'); 
fn_path     = fullfile(fn_path, 'Reconstruct_sortType(CC)_(11-Jun-2021)');
data_type   = upper(data_type);
fn_template = 'reconstruct_%s_(11-Jun-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';       

switch data_type
    case 'SU'
        unit_list = [10, 25, 50, 103]
    case 'MUA'
        unit_list = [10, 25, 50, 103];      %[10, 25, 50, 100, 150, 241];    
    otherwise
        error('--> Unrecognized DATA_TYPE!');
end

n_unit_list = length(unit_list);
drr = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
drr_idx = drr.sortby(1:n_drr);




%% Start Looping over all unit_list
for un = 1:n_unit_list
    % Loading first file to get preliminary metadata
    %       splits: [1×1 struct]
    %      spec_st: [1×1 struct]
    %      stim_st: [1×1 struct]
    %     tbl_data: [241×20 table]
    fn_n = sprintf(fn_template, data_type, unit_list(un));
    aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
    warning off
    dummy = load( fullfile(fn_path, fn_n), 'splits', 'spec_st', 'stim_st', 'tbl_data' );
    warning on

    % Extract data: 
    splits   = dummy.splits;
    spec_st  = dummy.spec_st;
    stim_st  = dummy.stim_st;
    tbl_data = dummy.tbl_data;
    n_splits = splits.n_splits;
    n_bands  = spec_st.n_bands;
    n_time   = length(splits.idx);
    % n_lags: # of samples along the x-axis (time) for each patch to use for comparison

    % Sampling frequency along the time axis
    binwidth    = spec_st.binwidth;     % (ms)
    fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
    f           = spec_st.f;            % (Hz)
    win_size_ms = spec_st.win_size_ms;  % (ms)
    t           = spec_st.t;            % (sec)


    %%
    dummy = nan(n_drr, n_splits);
    tbl.CC = array2table( dummy );
    tbl.CC.Properties.VariableNames = arrayfun(@(N) sprintf('sp%d', N), 1:n_splits, 'uni', 0);
    tbl.CC.Properties.RowNames = drr.labels(drr.ordered);

    % Fast duplicate, MATLAB style
    tbl.CC2 = tbl.CC;
    tbl.CC3 = tbl.CC;


    %% CC vs frequencies
    dummy = cell(n_drr, n_splits);
    tbl.CCf = array2table( dummy );
    tbl.CCf.Properties.VariableNames = arrayfun(@(N) sprintf('sp%d', N), 1:n_splits, 'uni', 0);
    tbl.CCf.Properties.RowNames = drr.labels(drr.ordered);

    % Fast duplicate, MATLAB style
    tbl.CCf2 = tbl.CCf;
    tbl.CCf3 = tbl.CCf;


    % % Pearson correlation for average of ALL the split
    % CC = nan(n_drr, n_splits);          % Sdry-vs-Sest
    % CC2 = nan(n_drr, n_splits);         % Sdrr-vs-Sest
    % CC3 = nan(n_drr, n_splits);         % Sdry-vs-Sdrr

    % Instant Pearson correlation as a function o DRR
    CCt = nan(n_time, n_drr);    	% Sdry vs Sest
    PPt = nan(n_time, n_drr);

    CCt2 = nan(n_time, n_drr);     	% Sdrr vs Sest
    PPt2 = nan(n_time, n_drr);

    CCt3 = nan(n_time, n_drr);   	% Sdry vs Sdrr ;only STIMULI
    PPt3 = nan(n_time, n_drr);



    % Concatenate all the estimations into one big matrix 
    Sest = nan(n_bands, n_time, n_drr);
    assert( 0 == sum(sum(isnan(spec_st.Sft{1}))), '--> There are NaNs in the DRY spectrograms!' );


    % Load the file
    fn_n = sprintf(fn_template, data_type, unit_list(un));
    aux.cprintf('Comments', '\n-> Loading <%s>...\n', fn_n);
    warning off
    data = load(fullfile(fn_path, fn_n));
    warning on
    obj_list = data.obj_list;
    spec_st  = data.spec_st;
    splits   = data.splits;

    % # of neurons
    assert( unit_list(un) ==  obj_list{1}.n_neurons, 'ERROR: something is wrong!' );
    n_neurons = obj_list{1}.n_neurons;

    if verbose
        fprintf('\n--> data_type  : %s\n', data_type);
        fprintf('--> n_neurons  : %d\n', n_neurons);
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_splits   : %d\n', n_splits);
    end


    debug_flag = 0;
    clc

    % Loop over the splits (exclusive intervals)
    for sp = 1:n_splits 
        idx_sp = sp == splits.idx;

        % Cut out the testing chunk of the spectrogram
        Sdry_ = spec_st.Sft{drr.dry}(:, idx_sp);    

        for k = 1:n_drr    
            rv = drr.ordered(k);

            % Cut out the testing chunk of the spectrogram
            Sdrr_ = spec_st.Sft{rv}(:, idx_sp);           
            Sest_ = obj_list{k,sp}.X_est;
            assert(0 == nnz(isnan(Sest_)))        

            % Concatenate all reconstructions\estimations into one big 3D
            Sest(:, idx_sp, k) = Sest_;

            % AVERAGED Pearson correlation coefficient:
            % DRY-vs-EST
            gof = goodness(Sdry_, Sest_);  
            tbl.CC{k,sp} = gof.CC;
            tbl.CCf{k,sp} = {gof.CCf};

            % DRR-vs-EST
            gof2 = goodness(Sdrr_, Sest_);   
            tbl.CC2{k,sp} = gof2.CC;
            tbl.CCf2{k,sp} = {gof2.CCf};

            % DRR-vs-EST
            gof3 = goodness(Sdry_, Sdrr_);   
            tbl.CC3{k,sp} = gof3.CC;
            tbl.CCf3{k,sp} = {gof3.CCf};

            % INSTANTANEOUS Pearson correlation coefficient:
            %CCt(idx_sp,k) = mean(Adry .* Aest)';
            [CCt(idx_sp,k),  PPt(idx_sp,k)]  = corrcoef_array(Sdry_, Sest_);
            [CCt2(idx_sp,k), PPt2(idx_sp,k)] = corrcoef_array(Sdrr_, Sest_);
            [CCt3(idx_sp,k), PPt3(idx_sp,k)] = corrcoef_array(Sdry_, Sdrr_);

            % DEBUG:
            if debug_flag
                figure(99);
                clf;
                tidx = t(idx_sp);
                plot(tidx, [CCt(idx_sp,k), CCt2(idx_sp,k)] );
                legend('CC(Dry-est)', sprintf('CC(%s-est)', drr.labels{rv}));
                xlabel('Time (sec)');
                ylabel('CC');
                title(sprintf('%s %d', data_type, unit_list(un)));
            end
        end
    end


    assert( 0 == sum(isnan(Sest(:))), ...
        '--> There are NaNs in one (or more) of the reconstructed spectrograms!' );

    info = ['- All DRR dimensions are already SORTED!\n',...
        '- Sest: (n_bands, n_time, n_drr)\n',...
        '\n',...
        '- CC : Sdry-vs-Sest, averaged over time\n',...
        '- CC2: Sdrr-vs-Sest, averaged over time\n',...
        '- CC3: Sdry-vs-Sdrr, averaged over time\n',...
        '- CCt : Sdry-vs-Sest, as a function of time\n',...
        '- CCt2: Sdrr-vs-Sest, as a function of time\n',...
        '- CCt3: Sdry-vs-Sdrr, as a function of time\n'...
    ];



    %% Save the analysis results
    % %{
    fprintf('\n- About to save the analysis results...\n');
    FN_PATH = load.path_to_data('Analysis');
    fn_name = sprintf('analyzed_%s_CC_(%s)_unit_list(%d).mat', data_type, date, unit_list(un));
    fn_fullfile = fullfile( FN_PATH, fn_name );

    % Save the results for that 
    save(fn_fullfile, ... 
        ... '-v7.3', ...
        'info',...              general info about the data
        'splits', ...           saves the chunks\intervals\speakers
        'stim_st', ...          stimulus data
        'spec_st', ...          spectrogram's structue with all relevant data
        'tbl_data', ...         a table with all neurons in the data set            
        'Sest',...              reconstructed\estimated spectrograms for all DRRs             
        'tbl', ...'CC', 'CCt', 'PPt',...
        ...'CC2', 'CCt2', 'PPt2', ...
        ...'CC3', 'CCt3', 'PPt3' ...
        'CCt', 'PPt', 'CCt2', 'PPt2', 'CCt3', 'PPt3' ...
    ); 

    fprintf('--> File <%s> SAVED!\n', fn_fullfile);
    %}
end













