%
% test_load_raw.m
%

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;

fignum = 11;

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
fn.session     =  'OBZ_E16_002';                    % 'OBZ_C74_029';
fn.file        = [fn.session, '-2-2-Spch'];     % '-2-2-Spch'
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
channel   = 6;    %1+fix(spikechan/6);     % channels number
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
mua = zeros(rawdata.smp_duration/200, rawdata.n_trials);
lfp = zeros(rawdata.smp_duration/200, rawdata.n_trials);

for in = 1:rawdata.n_inner
    for tr = 1:rawdata.n_trials_per_stimulus
        % Load the desired RAW measurement
        x = load.raw_extract_meas(X, rawdata, channel, tr, in, outer);  
        
        % aMUA
        ii = tr + (in-1)*rawdata.n_trials_per_stimulus;
        [mua(:,ii), Fs_, t_] = aMUA(x, rawdata.Fs);
        
        % LFP
        lfp(:,ii) = LFP(x, rawdata.Fs);
    end
end

mua = (mua - mean(mua))./std(mua);  % Z-score


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














