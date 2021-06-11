% plot_ONE_example_for_DRRs.m


clc
fignum = 11;
verbose = 1;

setup_environment('../../');
  
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);




%% Load a list of 
fn_unit_list = fullfile( load.path_to_data('Analysis'), 'idx_MUA_good_sorted_unit_thr(0-7).mat');
dummy = load(fn_unit_list);
sorted_unit_list = dummy.unit_list;



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'SU';       % {'SU', MUA'}
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
        unit_list = 103; %[10, 25, 50, 103];

    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';  
        unit_list = 241; %[10, 25, 50, 103, 241];

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

fn.load.path = load.path_to_data('data');
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data = load(fn.load.fullfile);

spec_st   = data.spec_st;
stim_st   = data.stim_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));

n_units   = height(tbl_data);     % total available units
len_unit_list = length(unit_list);
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
duration_sec   = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

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




%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( spec_st );



%% Splits the overall stimulus into chunks according to the speakers
binwidth       = spec_st.binwidth;     % (ms)

[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);



%%
switch upper(data_type)
    case 'SU'
        figure(0+fignum);  
    case 'MUA'
        figure(5+fignum);        
end

fontsize = 24;
fontsize_big = 32;
fontsize_bigger = 32;
markersize = 30;


un = 20 %11;         % selected unit to plot
spk = 12;        % speaker 1-12

Y = data.H(:, drr.ordered, un);
Y = Y./max(Y);
t = data.spec_st.t;             % (ms)
ax = zeros(1, n_drr);
idx = split_time_idx(spk,1):split_time_idx(spk,2);
ymax = -inf;
ymin = inf;


for k = 1:n_drr
    ax(k) = subplot(n_drr, 1, n_drr-k+1);
    
    %plot(t(idx) - t(idx(1)), Y(idx, k), 'Color', aux.rpalette(k));    
    switch upper(data_type)
    case 'SU'
        barh = bar(t(idx) - t(idx(1)), Y(idx, k), 5);
        set(barh, 'FaceColor', aux.rpalette(k));
        %set(barh, 'FaceAlpha', 0.8);    
    case 'MUA'
        ybias = 1.0;        
        plth = plot(t(idx) - t(idx(1)), 0.25*zca(Y(idx, k))+ybias);
        set(plth, 'Color', aux.rpalette(k));
        set(plth, 'LineWidth', 3);
    end

    set(ax(k), 'fontsize', fontsize);

    ymax = max( [ymax, ylim] );
    ymin = min( [ymin, ylim] );
    
    ylabel( drr.labels{drr.ordered(k)}, 'fontsize', fontsize_big);
    
    if 1<k 
        set(ax(k), 'XTickLabel', '');
    end
end


% linkaxes(ax);
title_str = aux.ctitle('Selected Response of a Female Speaker Uttering', ['"', tbl_metadata.txt{spk}(9:end-2), '"']);
title(ax(end), title_str, 'fontsize', fontsize_bigger);
xlabel(ax(1), 'Time (sec)', 'fontsize', fontsize_big);
axis(ax, 'tight')
set(ax, 'YTickLabel', '');































