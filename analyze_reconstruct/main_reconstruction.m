%
% main_reconstruction.m
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

 
 
 %% Load data
data_type = 'MUA';     % {'PSTH', 'MUA'}

if ~exist('H', 'var')
    aux.cprintf('g', '-> Loading NEW %s data...\n', data_type);
    %
    %            H: [18000󬝯56 double]
    %     H_labels: [356�6 double] (holds boolean labels, 1 means this
    %                              the measurement is available; 0 not.
    %      neurons: [1�437 double]
    %      spec_st: [1�1 struct]
    %      tbl_slc: [356�20 table]
    
    fn.path = '../.data/';
    
    switch upper(data_type)
        case 'PSTH'
            %fn.load = 'data_PSTH_bw(10)_fbands(30)_win(25)ms_(18-May-2020).mat';
            %fn.load = 'data_PSTH_bw(2)_fbands(30)_win(10)ms_(02-Jun-2020).mat';
            %fn.load = 'data_PSTH_bw(10)_fbands(30)_win(10)ms_(02-Jun-2020).mat';
            fn.load = 'data_PSTH_(08-Jun-2020)_bw(10)_fbands(30)_win(5)ms_spec(gammatone).mat';
        case 'MUA'
            %fn.load = 'data_MUA_bw(10)_fbands(30)_win(25)ms_(18-May-2020).mat'; 
            %fn.load = 'data_MUA_bw(2)_fbands(30)_win(5)ms_(02-Jun-2020).mat'; 
            %fn.load = 'data_MUA_bw(10)_fbands(30)_win(5)ms_(02-Jun-2020).mat'; 
            fn.load = 'data_MUA_(08-Jun-2020)_bw(10)_fbands(30)_win(5)ms_spec(gammatone).mat';
            
        otherwise
            error('-> Unidentified DATA_TYPE!');
    end
    load(fullfile(fn.path, fn.load));
    H0 = H;
    
else
    H = H0;      
    aux.cprintf('r', '-> *** Did NOT load new data! ***\n');
  
end
 
aux.cprintf('g', '-> Data file: %s!\n', fn.load);


binwidth    = spec_st.binwidth;     % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 

if verbose
    fprintf('--> binwidth   : %d\n', binwidth);
    fprintf('--> n_bands    : %d\n', n_bands);
    fprintf('--> win_size_ms: %d\n', win_size_ms);
end


%% Loads the CFs
%{
%   Name              Size            Bytes  Class     Attributes
% 
%   CF              437x6             20976  double              
%   Tpeak           437x6             20976  double              
%   best_rate       437x6             20976  double              
%   best_scale      437x6             20976  double              
%   cf_width        437x6             20976  double              
%   unit_rows       437x6             20976  double    
'DEBUG: change this; the data in here should be included in the usual loaded file'
load('..\.data\analyze_STRF-pars_MUA_bw(10)_fbands(30)_win(50)ms_(01-Jun-2020).mat');

switch upper(data_type)
    case 'MUA'
        % A list of invariant-CF unit indices (up to 1 semitone)
        % This can be calculated in **analyze_strf_CF_invariant.m**
        neu_cfinv = [7, 12, 21, 22, 24, 25, 31, 59, 67, 79, 100, 103, 154, ...
             163, 186, 187, 253, 255, 263];
    otherwise
        H_labels = ones(size(H0,3), 6);
        neu_cfinv = [];
end
%}
 

%% Selecting neurons to work with
% %{
'*** DEBUG: Selecting neurons to work with'

        % Get the valid measurements\columns
        invalid_neurons = arrayfun(@(I) 5 > nnz(~isnan(sum(H(:,:,I))) |...
            0 == sum(H(:,:,I))), 1:size(H,3));
        H = H(:,:,~invalid_neurons);


n_neurons = 25;
neurons = randperm(size(H,3), n_neurons);
% neurons = neu_cfinv;
% neurons = [neurons, randperm(size(H,3), 5)];     '^%^%%%$ DEBUG ^^%^#$ !!!'
% n_neurons = length(neurons);

H = H(:,:,neurons);
H_labels = H_labels(neurons,:);

% Remove empty entries
invalid_neurons = 0 == squeeze(sum(sum(H)))';
H = H(:,:,~invalid_neurons);
H_labels = H_labels(~invalid_neurons,:);
 
if verbose
    fprintf('--> n_neurons  : %d\n', n_neurons);
end
       
%}



%% Train-test split
aux.vprint(verbose, '\n--> Starting Train-test split\n');
train_drr = 3
test_drr  = 3
n_splits  = 24; 
test_grp_number = 1;
aux.vprint(verbose, '\t--> train group: %d\n', train_drr);
aux.vprint(verbose, '\t--> test group : %d\n', test_drr);

% Some of the measurements are missing
good_idx = logical( H_labels(:, train_drr) );

% Split for the TRAINING data
X1 = spec_st.Sft{train_drr};
y1 = squeeze( H(:, train_drr, good_idx) );

%         X1 = X1(7,:);  '######'
%     y1 = zca(y1);   '######'
%     y1 = max(y1, 0);

[X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
    'n_splits', n_splits, ...    
    'test_grp', test_grp_number ...
);


% Split for the TESTING data
X2 = spec_st.Sft{test_drr};
y2 = squeeze( H(:, test_drr, good_idx) );


%         X2 = X2(7,:);  '######'
%     y2 = zca(y2);     '######'
%     y2 = max(y2, 0);

[~, X_test, ~, y_test] = train_test_split(X2, y2, ...
    'n_splits', n_splits, ...    
    'test_grp', test_grp_number ...
);



%%
train.t = spec_st.t(splits.idx ~= test_grp_number);
train.f = spec_st.f;
train.CFs = tbl_impale.CF;
train.n_neurons = size(y1,2);
test.t = spec_st.t(splits.idx == test_grp_number);
test.f = spec_st.f;



%% RECONSTRUCTION
aux.vprint(verbose, '\n--> Start RECONSTRUCTION...\n');

jackknife_flag = 1

% Parameters for the ridge regression
% inv_type = [1, 0.5];
inv_type = [2, 0.005];

% Initialize the reconstruction object
obj = reconstruct(binwidth, 'inv_type', inv_type, 'f', spec_st.f); 

% Fit the model
obj.fit(X_train, y_train,...
    'jackknife_flag', jackknife_flag,...
    'n_splits', 23, ...
    'fignum', 11);


% Predict the spectrogram
X_est = obj.predict(y_test);

if verbose
    aux.cprintf('g', '--> RECONSTRUCTION SCORE: %g\n', goodness(X_test, X_est));
    fprintf('\n--> Finished RECONSTRUCTION\n');
end


%% Estimation errors
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



%% Save the data
%{
    'SAVE the analysis!'    
    fn.save = sprintf('../.data/reconstruct_MUA_bw(%g)_fbands(%d)_win(%g)ms_(%s).mat',...
        binwidth, n_bands, win_size_ms, date);
    obj_st = struct(obj);
    
    %save(fn.save);  % Save the whole workspace!
    save(fn.save, '-v7.3', 'X_test', 'X_est', 'H', 'spec_st', 'obj_st');
%}


%% Plot #1: the reconstruction
if isempty(fignum), return; end

% SUBPLOT #1:
figure(10+fignum);
clf

% time axis for the test spectrogram
binwidth = spec_st.binwidth;
dt = diff(spec_st.t(1:2));
t_test = [0:size(X_test,2)-1] * (1e-3*binwidth);
f_test = spec_st.f;

ax = subplot(2,1,1);
spec.plot_spectrogram(t_test, f_test, X_test, 10+fignum);
colorbar('off');
set(ax(1), 'XTickLabel', '');
xlabel('');
ylabel('');
title('$S_{test}(f,t)$')


ax(2) = subplot(2,1,2);
[~, surf_h] = spec.plot_spectrogram(t_test, f_test, X_est, 10+fignum);
colorbar('off');
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
obj.plot_reconstruction_filters(spec_st.f, 20+fignum);


%% Plot #4: Plot the conditional level densities (CLD)
figure(25+fignum);
clf;    

[px_y, C] = CLD_plot(X_test, X_est, 20, 25+fignum);








