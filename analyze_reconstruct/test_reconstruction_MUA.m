%
% main_reconstruction_MUA.m
%

clc
% close all
% clear all

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find(arrayfun(@(X) str2double(X), tbl_impale.SPK));
    spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;



%% Load MUA data
% !!! NOTE !!!
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
% Example of the <fn.load> file,
% 
%     struct with fields:
% 
%            H: [18000󬝯56 double]
%     H_labels: [356�6 double] (holds boolean labels, 1 means this
%                              the measurement is available; 0 not.
%      neurons: [1�437 double]
%      spec_st: [1�1 struct]
%      tbl_slc: [356�20 table]

% fn.load = fullfile('../.data', 'data_MUA_bw(2)_fbands(30)_win(30)ms_(07-May-2020).mat');
% fn.load = fullfile('../.data', 'data_MUA_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone)');
fn.load = fullfile('../.data', 'data_MUA_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)');
% fn.load = fullfile('../.data', 'data_MUA_(23-Jul-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone).mat');
aux.vprint(verbose, '--> [main_reconstruction_PSTH.m]: Loading <%s> ...\n', fn.load);
load(fn.load);

% Remove units with NaNs
valid_idx = logical( squeeze(prod(~isnan(sum(H(:,1:drr.n_drr,:),1)),2)) );
H = H(:,:,valid_idx);

        %'^%^%%%$ DEBUG ^^%^#$ !!!'
        %neurons = randperm(size(H,3), 5);
        %H = H(:,:,neurons);
        %H_labels = H_labels(neurons,:);
        
        % Remove empty entries
        %invalid_neurons = 0 == squeeze(sum(sum(H)))';
        %H = H(:,:,~invalid_neurons);
        %H_labels = H_labels(~invalid_neurons,:);

binwidth    = spec_st.binwidth;    % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 

%% Train-test split
aux.vprint(verbose, '\n--> Starting Train-test split\n');
train_drr = 3;
test_drr  = 4;
n_splits  = 12; 
test_split= 1;
dim_x     = 2;
dim_y     = 1;
aux.vprint(verbose, '\t--> train group: %d\n', train_drr);
aux.vprint(verbose, '\t--> test group : %d\n', test_drr);

% Some of the measurements are missing
%good_idx = logical( H_labels(:, train_drr) );

% Split for the TRAINING data
X1 = spec_st.Sft{train_drr};
%y1 = squeeze( H(:, train_drr, good_idx) );
y1 = squeeze( H(:, train_drr, :) );


%     '#$#$$#$ DEBUG -- TIME LAG ##$#$'
%     I = 1:size(y1,1)-10;
%     delay = 5;
%     X1 = X1(:,I);
%     y1 = y1(delay+I, :);
%     whos y1 X1

    
[X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
    'n_splits', n_splits, ...    
    'test_grp', test_split, ...
    'dim_x', dim_x, ...
    'dim_y', dim_y ...
);


% Split for the TESTING data
X2 = spec_st.Sft{test_drr};
%y2 = squeeze( H(:, test_drr, good_idx) );
y2 = squeeze( H(:, test_drr, :) );


%     '#$#$$#$ DEBUG -- TIME LAG ##$#$'
%     X2 = X2(:,I);
%     y2 = y2(delay+I, :);


[~, X_test, ~, y_test] = train_test_split(X2, y2, ...
    'n_splits', n_splits, ...    
    'test_grp', test_split, ...
    'dim_x', dim_x, ...
    'dim_y', dim_y ...
);


train.t = spec_st.t(splits.idx ~= test_split);
train.f = spec_st.f;
train.CFs = tbl_impale.CF;
train.n_neurons = size(y1,2);
test.t = spec_st.t(splits.idx == test_split);
test.f = spec_st.f;



% RECONSTRUCTION
aux.vprint(verbose, '\n--> Start RECONSTRUCTION...\n');

lags_ms = 30;
jackknife_flag = 1;

% Parameters for the ridge regression
% inv_type = [1, 0.5];
inv_type = [2, 0.005];

% Initialize the reconstruction object
obj = reconstruct_c(binwidth,...
    'f', spec_st.f,...
    'iscausal', 0,...   
    'algo_type', 'regression',...  {'regression', 'asd'}
    'lags_ms', lags_ms ); 

% Fit the model
obj.fit(X_train, y_train,...
    ...'algo_type', algo_type,... 
    'jk_flag', jackknife_flag,...
    'n_splits', n_splits, ...
    'fignum', []);



% Predict the spectrogram
X_est = obj.predict(y_test);

gof = goodness(X_test, X_est);

if verbose
    aux.cprintf('g', '--> CC: %g\n', gof.CC);
    fprintf('\n--> Finished RECONSTRUCTION\n');
end




% Save the data
%{
    'SAVE the analysis!'    
    fn.save = sprintf('../.data/reconstruct_MUA_bw(%g)_fbands(%d)_win(%g)ms_(%s).mat',...
        binwidth, n_bands, win_size_ms, date);
    obj_st = struct(obj);
    
    %save(fn.save);  % Save the whole workspace!
    save(fn.save, '-v7.3', 'X_test', 'X_est', 'H', 'spec_st', 'obj_st');
%}



% Estimation errors
% Error between target & estimated interval
cov_maxlag = 1000;  % (samples)
err.CC     = lagged_corrcoef(X_test, X_est, cov_maxlag);
err.MSE    = mse(X_test - X_est);
err.nMSE   = 100*err.MSE/mse(X_test);

% Error between DRY & estimated interval
err.CC2dry  = lagged_corrcoef(X_test, X_est, cov_maxlag);
err.MSE2dry = mse(X_test0 - X_est);
err.nMSE2dry= 100*err.MSE2dry/mse(X_test0);


if verbose
    fprintf('\n--> Estimation errors\n');
    fprintf('       test label: %s\n', drr.labels{test_drr});
    fprintf('     ------------------------\n');    
    fprintf('       CC        : %g\n', err.CC);
    fprintf('       MSE       : %g\n', err.MSE);
    fprintf('       nMSE      : %.2f %%\n', err.nMSE);
    fprintf('     ------------------------\n');    
    fprintf('       CCs 2 Dry : %g\n', err.CC2dry);
    fprintf('       MSE 2 Dry : %g\n', err.MSE2dry);
    fprintf('       nMSE 2 Dry: %.2f %%\n', err.nMSE2dry);
    fprintf('     ------------------------\n');        
end



%% Plot #1: the reconstruction
if isempty(fignum), return; end

% SUBPLOT #1:
figure(10+fignum);
clf

% time axis for the test spectrogram
dt = diff(spec_st.t(1:2));
t_test = [0:size(X_test,2)-1] * (1e-3*binwidth);
f_test = spec_st.f;

ax = subplot(2,1,1);
spec.plot_spectrogram(t_test, f_test, X_test, 10+fignum);
set(ax(1), 'XTickLabel', '');
xlabel('');
ylabel('');
title('$S_{test}(f,t)$')


ax(2) = subplot(2,1,2);
[~, surf_h] = spec.plot_spectrogram(t_test, f_test, X_est, 10+fignum);
title('');
% colorbar('off');
drawnow;
ax(1).Position(3) = ax(2).Position(3);
title('$S_{est}(f,t)$')


% Set the color axis to be the same for both spectrograms
cmin = 0;               % min([X_test(:); X_est(:)]);
cmax = 0.5*(max(X_est(:)) + max(X_test(:)));   % max([X_test(:); X_est(:)]);
caxis(ax(1), [cmin, cmax])
caxis(ax(2), [cmin, cmax])

linkaxes(ax);


%% Plot #2: count CF
% SUBPLOT #1:
figure(15+fignum);
clf;

% SUBPLOT #4: show the CF of the used neurons (on log10 scale)
[counts, centers] = hist(train.CFs, train.n_neurons);
bar(log10(centers), counts);   % log scale
ax = gca;
xlim(log10([2000, 18e3]));
x_ticks = get(gca, 'XTick');
x_new = arrayfun(@(X) num2str(1e-3*10.^X, '%.0f'), x_ticks, 'UniformOutput', 0);
ax.XTickLabel = x_new;
xlabel('CFs (kHz; $log_{10}$ scale)');
ylabel('Count');
title(sprintf('CFs of Used Neurons $(\\times%d)$', train.n_neurons));
ylim([0, 1.1*max(counts)]);


% Add ABCD labels
% aux.abc(ax, 48, 'northwestoutside');



%% Plot #3: plot the G filters
obj.plot_reconstruction_filters('f', 1e-3*spec_st.f, 'fignum', 20+fignum);


%% Plot #4: Plot the conditional level densities (CLD)
figure(25+fignum);
clf;    

[px_y, C] = CLD_plot(X_test, X_est, 20, 25+fignum);


%% Plot #5: correlations
%{
figure(30+fignum);
clf;

plot(1e-3*spec_st.f, sum(spec_st.Sft{2}/norm(spec_st.Sft{2}).*spec_st.Sft{3}/norm(spec_st.Sft{3}),2));
hold on
plot(1e-3*spec_st.f,  sum(X_test/norm(X_test).*X_est/norm(X_est),2), '--' )
hold off
xlabel('Frequency (kHz)');
%}

%% Plot #6: lagged autocorrelation between spectrogram frequency channels & PSTH 
% norm_cols = 0;
% ac = ac_spec_psth(spec_st.Sft{3}, squeeze(psth_mtx(:,3,:)), norm_cols, 35+fignum, 1e-3*spec_st.f);








