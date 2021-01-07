% plot_2_spectrograms.m
%
% Description:
% Use it for the slides to show one DRY and another reverberant (DRR = -8.2 dB)
% spectrograms.

clc
% close all
% clear all

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
addpath('../');

FigSetup;


%% Load data
% !!! NOTE !!!
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
% Load:
% 
% 'H', 'tbl_impale', 'spec_st', 'stim_st'
%
data_type    = 'MUA'    % {'SU', 'MUA'}
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
        %fn.load.file = 'data_SU_(26-Aug-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
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
data        = load(fn.load.fullfile);
spec_st     = data.spec_st;
binwidth    = spec_st.binwidth;     % (ms)
duration_ms = spec_st.duration_ms;
duration_sec= 1e-3*duration_ms;     % (sec)

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
aux.vprint(verbose, '-> data_type: %s\n', data_type);
        
 
%%
sp = 6;     % speaker # to plot
drr = get_DRR_list_and_indices;        
        
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, duration_sec);
     
% Time axis samples
n0 = split_time_idx(sp,1):split_time_idx(sp,2);
t  = (1e-3*binwidth)  * n0;
t = t - t(1);

f = spec_st.f;  % (Hz)

%% 
fontsize = 48;

% Spectrogram #1:
figh = figure(fignum);
clf;
X_dry = spec_st.Sft{drr.dry}(:, n0);    
spec.plot_spectrogram(t, 1e-3*f, X_dry, 'fignum', figh, 'fontsize', fontsize);
xlabel('Time (sec)');
ylabel('Frequency (kHz)');
% title(sprintf('%s', drr.labels{drr.dry}));
% set(gca, 'FontSize', ceil(0.75*fontsize));
ax = gca;
ax.CLim = [0, 90];
pos = get(gcf, 'Position');         % !!! FIGURE !!!


%
% Spectrogram #2:
figh = figure(5+fignum);
clf;
X_drr = spec_st.Sft{drr.ordered(end)}(:, n0);            
spec.plot_spectrogram(t, 1e-3*f, X_drr, 'fignum', figh, 'fontsize', fontsize);
xlabel('Time (sec)');
ylabel('Frequency (kHz)');
% title(sprintf('DRR: %s', drr.labels{drr.ordered(end)}));
% set(gca, 'FontSize', ceil(0.75*fontsize));
ax = gca;
ax.CLim = [0, 90];
set(gcf, 'Position', pos);     % !!! FIGURE !!!












