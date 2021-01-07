%
% test_find_best_subsets.m
%
% Description:
% Tests the reconstruction class.
%

clc
% close all
% clear all

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_impale', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;
n_drr = 5;  % N DRRs TO USE



%% Load data
% !!! NOTE !!!
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
% Load:
% 
% 'H', 'tbl_impale', 'spec_st', 'stim_st'
%
data_type    = 'MUA'    % {'SU', MUA'}
fn.load.path = '../.data';
switch upper(data_type)
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        
        %fn.load.file = 'data_SU_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
        fn.load.file = 'data_SU_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        
        %fn.load.file = 'data_MUA_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
        fn.load.file = 'data_MUA_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        %fn.load.file = 'data_MUA_(23-Jul-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone)';
        %fn.load.file = 'data_MUA_(03-Aug-2020)_bw(5)_fbands(60)_win(NaN)ms_spec(gammatone)';
        %fn.load.file = 'data_MUA_(03-Aug-2020)_bw(5)_fbands(30)_win(10)ms_spec(stft)';
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data    = load(fn.load.fullfile);
spec_st = data.spec_st;


aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
aux.vprint(verbose, '-> data_type: %s\n', data_type);



%% Get the valid measurements\columns
duration_sec = 36;  % (sec) stimulus duration to use

% Select all measurement with the desired duration
slc.duration = tbl_impale.duration_sec == duration_sec;

% Select all measurements with that contains all (5) sessions
% slc_valid_neurons = ~arrayfun(@(I) 5 > nnz(~isnan(sum(H(:,:,I))) | 0 == sum(H(:,:,I))), 1:size(H,3));
slc.valid_neurons = 5 <= sum( ~isnan( squeeze(sum(data.H,1)) ), 1)';
slc.valid_neurons = slc.valid_neurons(:);

% Indices of both boolean conditions
slc.valid_neuron_idx = slc.valid_neurons & slc.duration(data.neuron_list);


% # of units to use in this script
'**********************'
% [25 40 45 50 48]  %[20 32 10 15] %[2 30 31 32 35]  %[1 6 7] %[10, 15]
%  * [5 6 7 8 15 20]
% ** [67, 238, 161, 159, 228] 
% At first, include ALL neurons for the orthogonal analysis
units = 1:nnz(slc.valid_neuron_idx); 
'**********************'

% Choosing # of units
H_valid = data.H(:,:,slc.valid_neuron_idx);
H_unit = squeeze( H_valid(:,:,units) );


% A list of all "valid" units to use
slc.neuron_list = data.neuron_list(slc.valid_neuron_idx);
slc.neuron_list = slc.neuron_list(units);  % truncate the list



%%
iscausal    = 0;                    % use reconstructed causal filters? 
lags_ms     = 50;                 % (ms) maximum lags
binwidth    = spec_st.binwidth;     % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 
jk_flag     = 1;

% Computation parameters (Reconstuction object)
algo_type   = 1;    % 1: use convmtx, 2: use for-loop (Nima's lab) procedures
gpu_flag    = 0;    %0 < gpuDeviceCount;   % use the GPU, if available

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife: %d\n', jk_flag);
end



%%
% Option #2: 
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));

% Get the longest time interval
n_smp_split = max( diff(split_time_idx, [], 2) );



%% Train: orthogonalize the measurements and get the most effective set of units
train_drr   = 3     % drr.dry;
test_drr    = 3     % drr.dry;
test_speaker= 4;    % chunk\speaker 
n_units     = 5;    % Choose # of units to use for the reconstruction

aux.cprintf('string', '--> train_drr   : %d\t (%s)\n', train_drr, drr.labels{train_drr});
aux.cprintf('string', '--> test_drr    : %d\t (%s)\n', test_drr, drr.labels{test_drr});
aux.cprintf('string', '--> test_speaker: %d\t (%s)\n', test_speaker, tbl_metadata.fn{test_speaker});

% Get an optimal set
H_ = squeeze( H_unit(:, train_drr,:) );
[slc.optimal_sorted, P] = find_best_unit_set(H_);


% Use these units
%slc.optimal_sorted = slc.optimal_sorted(:)';
slc.optimal_units_used = slc.optimal_sorted(1:n_units, :);
slc.optimal_units_used = slc.optimal_units_used(:);

% Training responses
y1 = squeeze( H_unit(:, train_drr, slc.optimal_units_used) );

% Training spectrogram
X1 = spec_st.Sft{train_drr};

% Plot the projections 
if ~isempty(fignum)
    fprintf('--> slc.optimal_units_used: ');
    disp(slc.optimal_units_used(:)');

    figure(fignum);
    clf;
    plot(P, '.-');
    xlabel('Unit Number');
    ylabel('$Y^T\cdot Z$');
    title('Projections over Orthogonal Axes');
end



% Train the model using the "best" chosen units
[X_train, X_test0, y_train, ~, ~] = train_test_split(X1, y1, ...
    'split_time_idx', split_time_idx,...
    'test_grp', test_speaker ...
);

% Initialize the reconstruction object
obj = reconstruct_c(binwidth,...
    'f', spec_st.f,...
    'iscausal', iscausal, ...
    'lags_ms', lags_ms ); 

% Fit the model
obj.fit(X_train, y_train,...
    'jk_flag', jk_flag,...
    'n_splits', n_splits, ...
    'fignum', []);



% Test 
% Testing spectrogram
X2 = spec_st.Sft{test_drr};

% Testing responses
y2 = squeeze( H_unit(:, test_drr, slc.optimal_units_used) );

[~, X_test, ~, y_test, ~] = train_test_split(X2, y2, ...
    'split_time_idx', split_time_idx,...
    'test_grp', test_speaker );


% PREDICTION; Predict the spectrogram using the selected responses
obj.predict(y_test);    
gof = goodness(X_test, obj.X_est);



% Plot stuff
if verbose
    t = linspace(0, 1e-3*binwidth*size(X_test0,2), size(X_test0,2));    
    nolabels = true;
    if train_drr == test_drr
        ysub = 2;
    else
        ysub = 3;
    end
    
    figure(7 + fignum);
    clf;
    subplot(ysub,1,1);
    ax = spec.plot_spectrogram(t, 1e-3*obj.f, X_test0, 7+fignum, nolabels);
    set(ax(1), 'XTickLabel', '');
    ylabel( aux.ctitle(sprintf('$X_{train}$ (DRR: %d)', train_drr), '\\') );
    title(sprintf('Speaker %d, CC: %.3f, NMSE: %.3f (BW: %g)', test_speaker, gof.CC, gof.nmse, binwidth));
    
    subplot(ysub,1,2);
    ax(2) = spec.plot_spectrogram(t, 1e-3*obj.f, obj.X_est, 7+fignum, nolabels);
    set(ax(2), 'XTickLabel', '');
    ylabel( aux.ctitle(sprintf('$\\hat{X}_{test}$ (DRR: %d)', test_drr), 'Log Frequency (kHz)' ));    

    if train_drr ~= test_drr
        xlabel(ax(2), '');
        
        subplot(ysub,1,3);
        ax(3) = spec.plot_spectrogram(t, 1e-3*obj.f, X_test, 7+fignum, nolabels);
        ylabel( aux.ctitle(sprintf('$X_{test}$ (DRR: %d)', test_drr), '\\') );
    end
    xlabel('Time (sec)');
    linkaxes(ax);
    
    
    fprintf('--> GOF:\n');
    disp(gof);
    fprintf('\n');

    n_lags = obj.n_lags;
    ax = obj.plot_reconstruction_filters('fignum', 20+fignum);    
    set(ax, 'FontSize', 24);
    %drawnow;
end  











