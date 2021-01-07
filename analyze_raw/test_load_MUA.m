%
% test_load_MUA.m
%

clc

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup


%% Indices & labels of the stimuli
% Idx.      1        2       3        4        5       6
% DRR      Dry  9.4 dB  4.8 dB  -2.5 dB  -8.2 dB      --
% Revb.   100%     80%     80%      20%      20%    100%           
% Dist.   1.5m    1.5m    3.0m     1.5m     3.0m    3.0m
idx = get_DRR_list_and_indices;



%% Load neurons from excel sheet 
if ~exist('tbl_all', 'var') || isempty(tbl_all)
    data_st.path	= load.path_to_data; 
    data_st.fn      = 'LabNoteBook.xlsx';
    data_st.operator= 'OBZ';
    data_st.ID      = 'C74';
    data_st.measType= 'Spch';   % {'Spch', 'ns_Spch_fc1kHz', 'ns_Spch_fc4kHz'}
    data_st.session = [];     
    data_st.unit    = [];      
    data_st.measNum = [];     

    tbl_all = load.data_table(data_st);
end



%% STIMULI: load stimuli & measurement data
% ---------------------------------------------------------
neuron_stim = 50;  % 8, 15
% ---------------------------------------------------------

spectrogram_type = 'stft';      % {['matlab'], 'stft', 'multitaper'}
% spectrogram_type = 'multitaper';      % {['matlab'], 'stft', 'multitaper'}
f_scale     = 'log';	% {['lin'], 'log', 'erb'}
n_bands     = 30;   	% (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 2;        % (ms) binwidth of the resulted spectrogram 

% (ms) length of the temporal window; applies only for {'stft', 'multitaper'}
win_size_ms = 10;

[stim_st, spec_st] = load.stimulus_and_spectrogram(tbl_all, neuron_stim, ...
    'spectrogram_type', spectrogram_type, ...
    'binwidth', binwidth,...
    'n_bands', n_bands, ...
    'f_scale', f_scale, ... {'lin', 'log'}
    'win_size_ms', win_size_ms, ...    
    'nw', [], ...          (default: 1.4) only for spectrogram_type == MULTITAPER
    'fignum', [] ...
    );

% Duration of the stimulus, in mili-seconds
duration_ms = units.sec2ms( stim_st.info.Duration );

if verbose
    fprintf('\n--> STIMULUS...\n');
    fprintf('\t--> Sampling rate (fs): %g kHz\n', 1e-3*stim_st.fs);
    fprintf('\t--> Duration          : %g (sec)\n', 1e-3*duration_ms);
end




%% Impale files (MAT files)

% OdedAlienware
% path_root_mat = 'D:\DataBases\myMeas\Rabbits\OBZ_C74_022\';

% Oded's external HD
path_root_mat = 'C:\Users\barzelo\Google Drive\codeOnCloud\.data\';


%% Raw files (0.et files)
% OdedAlienware
% path_root_raw = 'D:\DataBases\myMeas\Rabbits\!RAW\';

% Oded's external HD
path_root_raw = 'D:\oded_backups\myDataBases\DataBases 03-17-2020\myMeas\Rabbits\!RAW\';

% Apollo drive
% path_root_raw = '\\apollo\ent_archive\Delgutte_Archive\obz\Rabbits\C74\!RAW\et_OBZ_C74_020\';


%% Set the files to load
fn.session     =  'OBZ_C74_021';                    % 'OBZ_C74_029';
fn.file        = [fn.session, '-1-11-Spch_dry'];     % '-2-2-Spch'
seps           = regexp(fn.file, '-');
fn.session_dir = fn.file(1:seps(1)-1);

% Impale
fn.path.Impale = [path_root_mat, fn.session_dir, filesep];
fn.Impale      = [fn.path.Impale, fn.file, '.mat'];
assert(~isempty(dir(fn.Impale)), sprintf('--> Can''t find %s', fn.Impale));

% Raw file
fn.path.raw = [path_root_raw, 'et_', fn.session_dir, filesep];
fn.raw      = [fn.path.raw, fn.file, '.0.et'];
assert(~isempty(dir(fn.raw)), sprintf('--> Can''t find %s', fn.raw));


%% LOADs the Impale structure & RAW data
% Impale structure
warning off
S_ = load(fn.Impale);
warning on

% Raw wav
[X, rawdata] = load.raw_spch_wave(S_, fn.raw);



%% Get the spikes times from the Impale structure
% For this demonstration, pick up one SPIKECHAN to work with
%
% Example: 
% 7/3/2019, OBZ_C74_020-1-2-Spch_dry
% spikechan: 1; n_spikes: 747; CF: 3810; Rating: B; SPL: 80 dB
channel   = 2;    %1+fix(spikechan/6);     % channels number
spikechan = 1 + (channel-1)*6;
trial     = 1;      % trial number to look at
inner     = 1;
outer     = 1;
duration_sec = S_.stimChans{1}.Source.numTokens;     % (sec) Stimulus duration

% Load the desired RAW measurement
x = load.raw_extract_meas(X, rawdata, channel, trial, inner, outer);  

% Get a specific trial
data = medit.get_Spch_tps(S_, spikechan);
tps_ = data.tps{inner}( data.trials{inner} == trial );


fprintf('\n Spikes times\n');
fprintf('--> spikechan   : %d\n', spikechan);
fprintf('--> channel     : %d\n', channel);
fprintf('--> idx_stim    : %d\n', inner);
fprintf('--> duration_sec: %g\n', duration_sec);
fprintf('--> numel(tps_) : %d\n', numel(tps_));



%% Plot the spike tarins
time_units = 'sec';     % {'sec', 'ms'}

figure(fignum);
clf;
t = linspace(0, duration_sec, rawdata.smp_duration)';   % (sec) time axis
if strcmpi('ms', time_units)
    t = 1e3*t;
end
plot(t, x);     % plot raw measurement
axis tight
hold on
% plot(1e-3*tps_, threshold_ch*ones(size(tps_)), 'v', 'markersize', 12);
spiketools.plot_dotraster(nan, data.t{inner}, data.ch{inner},...
    'spikechan', spikechan, ...
    'duration_sec', duration_sec,...
    'raster_type', 'lines',...
    'trial_to_plot', trial, ...
    'time_units', time_units );
xlabel(sprintf('Time (%s)', time_units));
hold off

% Plot spike threshold
threshold_V = S_.discrimSettings.settings(channel).threshold_V;
low_threshold = S_.discrimSettings.settings(channel).low_threshold;

% aux.hline(threshold_ch, 'LineWidth', 0.5);
aux.hline(threshold_V, 'LineWidth', 0.5);
aux.hline(low_threshold, 'LineWidth', 0.5, 'LineStyle', '-.');



%% aMUA
Fs = rawdata.Fs;        % (Hz)
Fs_new = 1/(1e-3*binwidth);    % (Hz)

% resampled stimulus
[Ny, Dy] = rat(Fs_new/stim_st.fs);
Y_ = resample(stim_st.Y, Ny, Dy);

% downsampling factors
[Nd, Dd] = rat(Fs_new/Fs);

mua = zeros(rawdata.smp_duration*(Nd/Dd), rawdata.n_trials);
lfp = zeros(rawdata.smp_duration*(Nd/Dd), rawdata.n_trials);


for in = 1:rawdata.n_inner
    for tr = 1:rawdata.n_trials_per_stimulus
        % Load the desired RAW measurement
        x = load.raw_extract_meas(X, rawdata, channel, tr, in, outer);  
        
        % aMUA
        ii = tr + (in-1)*rawdata.n_trials_per_stimulus;
        [mua(:,ii), Fs_, t_] = aMUA(x, Fs, Fs_new);
        
        % LFP
        lfp(:,ii) = LFP(x, Fs, Fs_new);
    end
end

% mua = (mua - mean(mua))./std(mua);  % Z-score


%% Plot
figure(5+fignum);
clf;

disp(500*1/Fs_);
[Sy, fy, ty] = spectrogram(mua(:,1), 500, 450, 60-1, Fs_,'yaxis');
imagesc(ty, fy, abs(Sy));
colorbar;
title('MUA');

figure(10+fignum);
clf;
[Slfp, fy, ty] = spectrogram(lfp(:,1), 500, 450, 60-1, Fs_,'yaxis');
imagesc(ty, fy, abs(Slfp));
colorbar;
title('LFP');



%%
% binwidth        = 1;   % (ms)
pstw_win_size_ms= 10;
syncchan        = 0;

%S_.info.duration_ms = duration_ms;
%S_.info.n_sync      = [1, 1, 1, 1, 1];
% pst_type = {'pstw', 30};    % {'pstw', WIN_SIZE_MS}    
% psth_mtx = calc_PSTHs(S_, binwidth, 'pst_type', pst_type);

data = medit.get_Spch_tps(S_, spikechan, syncchan);

psth_k = PSTW(data.t(1), data.ch(1),...
    duration_ms,...
    binwidth,...   % !! set the PSTH bins to equal that of the spectrogram
    pstw_win_size_ms, ...
    spikechan,...
    syncchan );

figure(11);
clf;
plot(1e-3*psth_k.bins, [psth_k.Hbin, psth_k.H]);

hold on
spiketools.plot_dotraster(nan, data.t{inner}, data.ch{inner},...
    'spikechan', spikechan, ...
    'duration_sec', duration_sec,...
    'raster_type', 'bars',...
    'trial_to_plot', trial, ...
    'time_units', time_units );
hold off


%%








