% plot_one_unit.m
%
% Plots one selected unit (single unit or multiunit activity) for the 
% presentation.

clc
fignum = 11;
verbose = 1;
%addpath('../');
FigSetup;



%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the mea
if ~exist('tbl_impale', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;
%n_drr = 5;  % N DRRs TO USE


%% META-DATA 
dummy = load('..\.data\stimulus\metadata_(36)_wav_(30-Jun-2020).mat');
tbl_metadata = dummy.tbl_metadata;



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
neuron_duration_list = tbl_impale.duration_sec == duration_sec;

% Select all measurements with that contains all (5) sessions
slc.valid_neuron_idx = 5 <= sum( ~isnan( squeeze(sum(data.H,1)) ), 1)';
slc.valid_neuron_idx = slc.valid_neuron_idx(:);

% Indices of both boolean conditions
slc.valid_neuron_idx = slc.valid_neuron_idx & neuron_duration_list(data.neuron_list);

% A list of all "valid" units to use
slc.valid_neurons = data.neuron_list(slc.valid_neuron_idx);
n_units = nnz(slc.valid_neuron_idx);

% Choosing # of units
H = data.H(:,:,slc.valid_neuron_idx);

% % For SUs, use a LPF FIR window
% H = data.H;
% if strcmpi('SU', data_type)
%     aux.cprintf('UnterminatedStrings', '--> Using a LPF GAUSSIAN window over the SUs!!!:\n');
%     win_lpf = gausswin(5);
%     H = filtfilt(win_lpf, 1, H);
% end


%% Load the STRFs data of SUs
% whos
%   Name               Size                  Bytes  Class     Attributes
% 
%   H_valid         7200x6x103            35596800  double              
%   scores             1x1                  149976  struct              
%   slc                1x1                    1326  struct              
%   spec_st            1x1                15039301  struct              
%   splits             1x1                   58592  struct              
%   tbl_impale       437x20                 640191  table               
%   tbl_strf         103x5                  373423  table               
fn.path = '..\.data\STRFs\';
fn.name = 'STRF_SU_(25-Aug-2020)_units(103)_bw(5)ms_algo(regression)_fbands(30)_lags(30)ms_cau(1)_trainDRR(3).mat';
% fn.name = 'STRF_SU_(27-Aug-2020)_units(103)_bw(1)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(3)';
fn.file2load = fullfile(fn.path, fn.name);
dummy = load(fn.file2load);
tbl_strf = dummy.tbl_strf;


%% Prepare the data to show
% Choose the unit to use
%
% Row #  BF(STRF)
%     6  7255 Hz
%   101  4254 Hz
%    45  5633 Hz
%    52  1264 Hz *CC non-decreasing
%    58  5870 Hz *
%    64  4009 Hz *
%
N          = 5% 53                     % unit # to use  
un         = slc.valid_neurons(N);  % measurement #
test_speaker= 6;                    % test speaker (1 to 12)
binwidth    = spec_st.binwidth;     % (ms)
duration_ms = spec_st.duration_ms;  % (ms)

% What is the speaker's sex?
if contains(tbl_metadata.fn{test_speaker}, '_M')
    test_sex = 'Male';
else
    test_sex = 'Female';
end

% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*duration_ms);


% BEST FREQUENCT
% assert( numel(slc.valid_neurons) == nnz(slc.valid_neurons == tbl_strf.unit) );
% bf = tbl_strf.bf(un == tbl_strf.unit);
% assert(~isempty(bf));
%     % [~, freq_idx] = max(spec_st.Sft{drr.dry} * H(:,drr.dry,N));
%     % bf = spec_st.f(freq_idx);
          bf = data.tbl_impale.CF(un);    % (Hz)

%         '^%^%#^$(*$(#$*(#%$*(%$'
%         bf = 6.2760e+03

fprintf('--> unit    : %d\n', un);
fprintf('--> BF(STRF): %g Hz\n', bf);


%
% ADD a plot of the broadband envelope
% - Load the DRY stimulus
dummy = load(fullfile('../.data/stimulus', 'data_stimuli_duration(36_and_40)sec.mat'));
stim_st = dummy.stim_list{spec_st.duration_ms == duration_ms};
fs = 1/(1e-3*binwidth);     % (Hz)
yenv = calc_stimulus_envelope(stim_st.Y(:,drr.dry), bf, fs, stim_st.fs);
y_idx = split_time_idx(test_speaker,:);
yenv_ = yenv(y_idx(1):y_idx(2));    % cut out the desired test speaker segment

%
figure(fignum);
clf;

fontsize = 24;
bar_width = 1.5;

nt    = diff(y_idx)+1;
t     = linspace(0, 1e-3*binwidth*nt, nt);
drr_inv = fliplr(drr.ordered);
ax = nan(1, drr.n_drr);
legendh = nan(1, drr.n_drr);
CC = nan(1, drr.n_drr);
for k = 1:drr.n_drr
    yk = H(y_idx(1):y_idx(2), drr_inv(k), N);
    ax(k) = subplot(6,1,k);
    color_k = aux.rpalette(sprintf('new%02d', drr.n_drr-k+1));
    
    if strcmpi('SU', data_type)
        bar(t, yk/max(yk), bar_width,...
            'FaceColor', color_k,...
            'EdgeColor', 'none');
    else
        plot(t, (yk-min(yk))/(max(yk)-min(yk)), 'Color', color_k);
    end
    
    ylabel(drr.labels{drr_inv(k)});
    CCk = corrcoef(yenv_, yk);
    CC(k) = CCk(1,2);
    legendh(k) = legend(ax(k), sprintf('CC: %.2f', CC(k)), 'Location', 'northeast');

end
% set(legendh, 'FontSize', ceil(0.7*fontsize));
set(ax(1:5), 'XTickLabels', '');
linkaxes(ax(1:5), 'y');
ylim([0.2, 1.0]);

% ### SUBPLOT 6 ###
ax(6) = subplot(6,1,6);
plot(t, yenv_/max(yenv_), '-k');
xlabel('Time (sec)');
ylabel(aux.ctitle('Dry','Envelope'));
linkaxes(ax, 'x');
axis tight
set(ax, 'FontSize', fontsize);
legendh(end+1) = legend(ax(6), '$y_{env}$', 'Location', 'northeast'); 
set(legendh, 'FontSize', floor(fontsize));

xlim([0,3.20]);


% set(ax, 'YTick', [0, 1]);
set(ax, 'YTickLabels', '');
title(ax(1), sprintf('%s Speaker (%s, $BF_{strf}$: %g Hz)', test_sex, data_type, round(bf)));

aux.abc(ax, 'fontsize', 34);


% ### FIGURE ###
% get(gcf, 'Position');
set(gcf, 'Position', [64, 2, 1366, 994]);


%% plot of CC vs. DRR
figure(5+fignum);
clf;

fontsize = 40;
markersize = 40;

CC_ = CC(end:-1:1);
plot(1:drr.n_drr, CC_, ':k', 'MarkerSize', markersize);
hold on
plth = nan(1, drr.n_drr);
for k = 1:drr.n_drr
    plth(k) = plot(k, CC_(k), 'o');
    set(plth(k), 'MarkerSize', markersize ,...
        'MarkerEdgeColor', aux.rpalette(k),...
        'MarkerFaceColor', aux.rpalette(k));
end
hold off
xlabel('DRR');
if strcmpi('SU',data_type)
    ylabel('CC(ENVELOPE-to-PSTH)');
elseif strcmpi('MUA',data_type)
    ylabel('CC(ENVELOPE-to-MUA)');
else
    error('');
end
set(gca, 'XTick', 1:drr.n_drr);
set(gca, 'XTickLabel', {drr.labels{drr.ordered}});
set(gca, 'FontSize', fontsize);

xlim([0.5, 5.5]);
ylim([0.95*min(ylim), 1.05*max(ylim)]);





