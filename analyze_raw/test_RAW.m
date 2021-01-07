%
% test_RAW.m
%
%
% Description:
% This file loads Impale structure and loads its raw data. 
% * The Impale files are MERGED files (S.info.merged == 1) that were created 
%    using StimViewerGUI. 
% * Raw data includes ALL recorded channels.
% 

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;

fignum = 11;


%% ROOT path
% OdedAlienware
% path_root_mat = 'D:\DataBases\myMeas\Rabbits\OBZ_C74_022\';

% Oded's external HD
path_root_mat = 'C:\Users\barzelo\Google Drive\codeOnCloud\.data\';


%% RAW WAVE files (0.et files)
% OdedAlienware
% path_root_raw = 'D:\DataBases\myMeas\Rabbits\!RAW\';

% Oded's external HD
path_root_raw = 'D:\oded_backups\myDataBases\DataBases 03-17-2020\myMeas\Rabbits\!RAW';

% Apollo drive
% path_root_raw = '\\apollo\ent_archive\Delgutte_Archive\obz\Rabbits\C74\!RAW\et_OBZ_C74_020\';


%% Impale + Path; 
% set the files to load
clear fn
fn.session     =  'OBZ_C74_019';                        % 'OBZ_E16_002'; 'OBZ_C74_029';
fn.file        = [fn.session, '-2-10-Spch_MERGED'];   	% '-2-2-Spch'  ; '-2-2-Spch'
seps           = regexp(fn.file, '-');
fn.session_dir = fn.file(1:seps(1)-1);

% Impale
fn.path.Impale = [path_root_mat, fn.session_dir, filesep];
fn.Impale      = [fn.path.Impale, fn.file, '.mat'];
assert(~isempty(dir(fn.Impale)), sprintf('--> Can''t find %s', fn.Impale));

% Impale structure
warning off
S_ = load(fn.Impale);
warning on


%% Load the RAW file
fn_raw_tbl = [path_root_raw, filesep, S_.info.rawDataDir, filesep, S_.info.rawFileName];
assert(~isempty(dir(fn_raw_tbl)), '--> ERROR: can''t find this file!!');
raw_st = load(fn_raw_tbl, '-mat');


%%
inner     = 3;
outer     = 1;
channel   = 1;    %1+fix(spikechan/6);     % channels number
spikechan = 1 + (channel-1)*6;
trial     = [];      % trial number to look at

aux.cprintf('blue', '\n--> File name: %s\n', fn.file);
aux.cprintf('blue', '--> inner: %d, outer: %d, channel: %d\n', inner, outer, channel);


% Get the table of all session measurements and add row number
T = raw_st.tbl;
T.row = [1:height(T)]';
T = [T(:,end), T(:,1:end-1)];

idx_tbl = (T.outer == outer) &...
    (T.inner == inner) &...
    (T.channel==channel);
assert(1<=nnz(idx_tbl), '--> [RAW waveform]: Did not find such measurement!');

n_trials = nnz(idx_tbl);

len_x = length(raw_st.tbl.x{1});    % (samples)

duration_ms = raw_st.rawdata.t_duration;
duration_sec = 1e-3*duration_ms;
t = linspace(0, duration_sec, length(T.x{1}))';   % (sec) time axis

% Extract the raw waveform vectors into one matrix
% X = [T.x{:}];  % this won't work if the measurement doesn't have all
%                % trials
X = nan(len_x, size(T,1));
for k = 1:size(T,1)
    % remove incomplete raw waveforms
    %if len(T.x{k}) ~= len_x
    xk = T.x{k};
    X(1:length(xk),k) = xk;
end



figure(fignum);
clf;
y_jumps = 1.5*max(max(abs(X(:,idx_tbl))));
plot(t, X(:,idx_tbl) + y_jumps*ones(len_x,n_trials)*diag([0:n_trials-1]));
set(gca, 'YTick', y_jumps*[0:n_trials-1]);
set(gca, 'YTickLabel', num2cell(1:n_trials));
ylabel('Trial');
xlabel(sprintf('Time (%s)', 'sec'));

axis tight
dist = unique(T.Dist(idx_tbl));     % (m)
if any(contains(T.Properties.VariableNames, 'Reverb'))
    reverb = unique(T.Reverb(idx_tbl));
elseif any(contains(T.Properties.VariableNames, 'Batch'))
    dist   = 1.5;   % (m) 
    reverb = unique(T.Batch(idx_tbl));
end
aux.ctitle(aux.mName2latex(fn.file), sprintf('Dist: %gm, Wall Coef: %g\\%%', dist, reverb));
% aux.ctitle(aux.mName2latex(fn.file), sprintf('Dist: %gm, Reverb: %g \\%%', dist, reverb));


%
% %{
trial = 1;      % trial number to look at
idx_tbl = (T.outer == outer) &...
    (T.inner == inner) &...
    (T.channel==channel) & ...
    (T.trial==trial);

x = X(:,idx_tbl);

figure(2+fignum);
clf;
plot(t, x);     % plot raw measurement
axis tight
hold on
% plot(1e-3*tps_, threshold_ch*ones(size(tps_)), 'v', 'markersize', 12);
spiketools.plot_dotraster(nan, S_.t{inner,outer}, S_.ch{inner,outer},...
    'spikechan', spikechan, ...
    'duration_sec', duration_sec,...
    'raster_type', 'lines',...
    'trial_to_plot', trial, ...
    'time_units', 'sec' );
xlabel(sprintf('Time (%s)', 'sec'));
ylabel('Amp.');
hold off
aux.ctitle(aux.mName2latex(fn.file), sprintf('Dist: %gm, Wall Coef: %g\\%%', dist, reverb));

% Plot spike threshold
threshold_V = S_.discrimSettings.settings(channel).threshold_V;
low_threshold = S_.discrimSettings.settings(channel).low_threshold;

% aux.hline(threshold_ch, 'LineWidth', 0.5);
aux.hline(threshold_V, 'LineWidth', 0.5);
aux.hline(low_threshold, 'LineWidth', 0.5, 'LineStyle', '-.');

%}


% spec.spectrogram(LFP(X(:,2),10e3), 500,...
%     'duration_ms', 40e3,...
%     'method', 'matlab',...
%     'lowfreq', 1, ...
%     'highfreq', 250,...
%     'binwidth', 2,...
%     'fignum', 10);

%%
%{
lfp = nan(18e3, size(X,2));
mua = nan(18e3, size(X,2));
color = nan(3, size(X,2));
for k = 1:height(T)
    mua(:,k) = MUA(X(:,k), raw_st.sr);
    lfp(:,k) = LFP(X(:,k), raw_st.sr);
    
    c = nan;
    if 1.5 == T.Dist(k) && 20 == T.Reverb(k)
        c = 4;
    elseif 1.5 == T.Dist(k) && 80 == T.Reverb(k)
        c = 2;
    elseif (1.5 == T.Dist(k) || 3.0 == T.Dist(k) ) && 100 == T.Reverb(k)
        c = 1;
    elseif 3.0 == T.Dist(k) && 20 == T.Reverb(k)
        c = 5;
    elseif 3.0 == T.Dist(k) && 80 == T.Reverb(k)
        c = 3;
    end
    color(:, k) = aux.rpalette(sprintf('new%02d', c));
end

% scatter(mean(mua), std(mua), 25, color', 'filled')
% scatter(skewness(mua), kurtosis(mua), 25, color', 'filled')

%}



