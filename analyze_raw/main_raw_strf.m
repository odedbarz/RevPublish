%
% main_RAW_STRF.m
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

% Whitening sphere
zca = @(x) (x-mean(x))./std(x);


%% Select a neuron to work with

% =================================
meas_row_number = 39;    
% =================================



%% Indices & labels of the stimuli
% Idx.      1        2       3        4        5       6
% DRR      Dry  9.4 dB  4.8 dB  -2.5 dB  -8.2 dB      --
% Revb.   100%     80%     80%      20%      20%    100%           
% Dist.   1.5m    1.5m    3.0m     1.5m     3.0m    3.0m
drr = get_DRR_list_and_indices;


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Add row numbers
    row = array2table( (1:height(tbl_impale))', 'VariableNames', {'row'} );
    tbl_impale = [row, tbl_impale];
end


%% Load the Impale structure
aux.vprint(verbose, '\n--> Load response PSTH_MTX data...\n');

% Select the desired stimulus (36 sec or 40 sec)
tbl_slc         = tbl_impale(meas_row_number, :);
duration_sec    = tbl_slc.duration_sec;     % (sec) stimulus duration
S               = load.response( tbl_slc, duration_sec );
S               = S{1};     % cell --> struct
duration_ms     = S.info.duration_ms;
spikechan       = tbl_slc.spikechan;
syncchan        = 0;



%% STIMULI: load stimuli & measurement data
%  clear stim_st

if  1       % ~exist('stim_st', 'var')
    fprintf('--> Load stimuli & spectrograms\n');
    spectrogram_type = 'stft';      % {['matlab'], 'stft', 'multitaper'}
    % spectrogram_type = 'multitaper';      % {['matlab'], 'stft', 'multitaper'}
    f_scale     = 'log';	% {['lin'], 'log', 'erb'}
    n_bands     = 30   	% (1x1) # of bins along the frequency domain of the spectrogram
    binwidth    = 10        % (ms) binwidth of the resulted spectrogram 
    win_size_ms = 25       % (ms) temporal window size over which to calc the spectrogram
    
    [stim_list, spec_list] = load.stimulus_and_spectrogram(tbl_impale, ...
        'spectrogram_type', spectrogram_type, ...
        'binwidth', binwidth,...
        'win_size_ms', win_size_ms, ...
        'n_bands', n_bands, ...
        'f_scale', f_scale, ... {'lin', 'log'}
        'nw', [], ...          (default: 1.4) only for spectrogram_type == MULTITAPER
        'fignum', [] ...
        );
    assert(spec_list{1}.binwidth == spec_list{2}.binwidth);
    
    % All available stimulus durations
    duration_all_sec = cellfun(@(X) X.info.Duration , stim_list, 'UniformOutput', 1);
    
    % Get the DURATION stimulus that corresponods to the desired selected neuron
    duration_sec= tbl_slc.duration_sec;
    
    % Make sure that all agrees
    assert(units.sec2ms( duration_sec ) == duration_ms);

    % Plot the spectrogram
    spec.plot_spectrogram(spec_list{1}.t, spec_list{1}.f, spec_list{1}.Sft{3});
    xlim([1,4]);
    
    idx_stimuli = find(duration_all_sec == duration_sec);
    stim_st     = stim_list{idx_stimuli};
    spec_st     = spec_list{idx_stimuli};       
end


if verbose
    fprintf('\n    STIMULUS info\n');
    fprintf('--> binwidth: %g ms\n', spec_st.binwidth);
    fprintf('--> Duration: %g sec\n', duration_sec);
end



%% Load the RAW structure file
fn.raw.fulfilename = [path_root_raw, filesep, S.info.rawDataDir, filesep, S.info.rawFileName];
assert(~isempty(dir(fn.raw.fulfilename)), '--> ERROR: can''t find this file!!\n');
raw_st = load(fn.raw.fulfilename, '-mat');

% Get the raw table. It contains all information and the raw waveforms
tbl_raw = raw_st.tbl;



%% Calculate MUA & LFP
clear mua
[mua.avg, mua.SE] = calc_raw_means(tbl_raw, raw_st.sr, binwidth, duration_sec, tbl_slc.spikechan, 'MUA');
H = mua.avg;

%clear lfp
%[lfp.avg, lfp.SE] = calc_raw_means(tbl_raw, raw_st.sr, tbl_slc.spikechan, 'LFP');
% H = lfp.avg;



%% Train\Test split
aux.vprint(verbose, '\n--> Starting Train-test split\n');
grp.train = 3;
grp.test  = 3;
n_splits  = 12; 
test_split= 5;
dim_x     = 2;
dim_y     = 1;

% Split for the TRAINING data
X1 = spec_st.Sft{grp.train};
y1 = H(:, grp.train);
[X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
    'n_splits', n_splits, ...    
    'test_grp', test_split, ...
    'dim_x', dim_x, ...
    'dim_y', dim_y ...
);

% Split for the TESTING data
X2 = spec_st.Sft{grp.test};
y2 = H(:, grp.test);
[~, X_test, ~, y_test] = train_test_split(X2, y2, ...
    'n_splits', n_splits, ...    
    'test_grp', test_split, ...
    'dim_x', dim_x, ...
    'dim_y', dim_y ...
);

train.idx = splits.idx ~= test_split;
train.t = spec_st.t( train.idx );
train.f = spec_st.f;
% train.CFs = tbl_impale.CF;
train.n_neurons = size(tbl_impale,1);

test.idx = splits.idx == test_split;
test.t = spec_st.t( test.idx  );
test.f = spec_st.f;

aux.vprint(verbose, '\t--> train DRR case: %s\n', drr.labels{grp.train});
aux.vprint(verbose, '\t--> test DRR case : %s\n', drr.labels{grp.test});


%% STRF
% (samples) the window to take for the xcorr (maxlags)
n_win = 140; %2*39+1;     

% Ridge regression constant(s)
tol = 0.005;

% calculate the strf(s)
[strf, strf_st] = strfpkg.strf(X_train, y_train, n_win,...
    'tol', tol,...
    'faxis', spec_st.f, ...
    'binwidth', spec_st.binwidth, ...     
    'fignum', []);


figure(10+fignum);     
%fignum = fignum + 2;
clf;
strfpkg.plot_strf(strf, spec_st.f, binwidth);


% Estimate the response
r_est = strfpkg.reconstruct_response( X_test, strf,...
    'avg', strf_st.avg.psth, ...
    'normalize', max(y2));

%
figure(15+fignum);     
plot(test.t, zca([y2(test.idx), r_est]));
ylabel('Response (whitened)');
xlabel('Time (sec)');
CC = corrcoef(y2(test.idx), r_est);
title(sprintf('CC: %.3f', CC(1,2)));
legend('$r_{test}$', '$\hat{r}$');
set(gca, 'FontSize', 20);




return;
'#######################'

%% STRF via optimized ridge regression
nks    = size(strf);       % number of filter pixels along [cols, rows]
nk     = prod(nks);        % total number of filter coeffs
nsamps = nnz(train.idx);

x_ = (X_train - mean(X_train,2))./std(X_train,[],2);
x = genstimhist(nks(2), x_);
y = zca( y_train );

% Sufficient statistics (old way of doing it, not used for ASD)
dd.xx = x'*x;   % stimulus auto-covariance
dd.xy = (x'*y); % stimulus-response cross-covariance
dd.yy = y'*y;   % marginal response variance
dd.nx = nk;     % number of dimensions in stimulus
dd.ny = nsamps;  % total number of samples

% Run ridge regression using fixed-point update of hyperparameters
% maxiter = 100;
lam0 = 100;
[kridge, hprs] = autoRidgeRegress_fp(dd, lam0);
fprintf('--> tol: %g\n', hprs.alpha*hprs.nsevar);

strf_ridge = reshape(kridge, nks);

% Estimate the response
r_est_ridge = strfpkg.reconstruct_response( X_test, strf_ridge,...
    'avg', strf_st.avg.psth, ...
    'normalize', max(y2));


% Plot
figure(20+fignum);     
clf;
strfpkg.plot_strf(strf_ridge, spec_st.f, binwidth);

figure(25+fignum);     
clf;
plot(test.t, zca([y2(test.idx), r_est_ridge]));
ylabel('Z-score');
CC = corrcoef(y2(test.idx), r_est_ridge);
title(sprintf('CC: %.3f', CC(1,2)));



%%

fprintf('\n\n...Running ASD_2D...\n');

% minlens = [2; 2];  % minimum length scale along each dimension
minlens = [1; 1];  % minimum length scale along each dimension

tic; 
[kasd,asdstats] = fastASD(x, y, nks, minlens);
toc;

strf_asd = reshape(kasd, nks);

% Estimate the response
r_est_asd = strfpkg.reconstruct_response( X_test, strf_asd,...
    'avg', strf_st.avg.psth, ...
    'normalize', max(y2));

% Plot
figure(30+fignum);     
clf;
strfpkg.plot_strf(strf_asd, spec_st.f, binwidth);

figure(35+fignum);     
clf;
plot(test.t, zca([y2(test.idx), r_est_asd]));
ylabel('Z-score');
CC = corrcoef(y2(test.idx), r_est_asd);
title(sprintf('CC: %.3f', CC(1,2)));





%% ISI
%{
% n_time   : # of time samples
% n_drr    : # of stimuli types (DRRs, like in the labels)
% n_neurons: # of measured neurons
[n_time, n_drr, n_neurons] = size(psth_mtx); 
[rates, sem] = calc_Rates(S_list, 0, duration_ms);                  % Rates & SEM
isi = calc_ISI(S_list, binwidth, 'n_bins', ceil(1000/binwidth));    % ISIs
%}

% whos isi rates psth_mtx   % DEBUG
aux.vprint(verbose, '--> Finished\n');















