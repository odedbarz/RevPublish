%
% main_reconstruction_SU.m
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
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find(arrayfun(@(X) str2double(X), tbl_impale.SPK));
    spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;


%% Load SUs data
% !!! NOTE !!!
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
% Example of the <fn.load> file,
% 
%     struct with fields:
% 
%            H: [18000×6×356 double]
%     H_labels: [356×6 double] (holds boolean labels, 1 means this
%                              the measurement is available; 0 not.
%      neurons: [1×437 double]
%      spec_st: [1×1 struct]
%      tbl_slc: [356×20 table]
fn.load = fullfile('../.data', 'data_SU_bw(2)_fbands(30)_win(20)ms_(07-May-2020).mat');
aux.vprint(verbose, '--> [main_reconstruction_SU.m]: Loading <%s> ...\n', fn.load);
load(fn.load);


%         '^%^%%%$ DEBUG ^^%^#$ !!!'
%         neurons = randperm(size(H,3), 10);
%         H = H(:,:,neurons);
% 
%         % Remove empty entries
%         invalid_neurons = 0 == squeeze(sum(sum(H)))';
%         H = H(:,:,~invalid_neurons);

binwidth    = spec_st.binwidth;    % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 


%% Train-test split
aux.vprint(verbose, '\n--> Starting Train-test split\n');
jackknife_flag= 1;
train_drr   = drr.dry;
% test_drr  = 3;
n_splits    = 20; 
test_grp_number =  5;  % set this group as a testing group out of the whole N_SPLITS groups
% n_drr       = length(drr.labels);
[n_smp, n_drr, n_neurons]= size(H);
n_smp_split = ceil(n_smp/n_splits);     % # of samples in each split


aux.vprint(verbose, '\t--> train group: %d\n', train_drr);
% aux.vprint(verbose, '\t--> test group : %d\n', test_drr);


%% Split the TRAINING set
X1 = spec_st.Sft{train_drr};
y1 = squeeze( H(:, train_drr, :) );
[X_train, X_test0, y_train, y_test0, splits] = train_test_split(X1, y1, ...
    'n_splits', n_splits, ...    
    'test_grp', test_grp_number ...
);


%% Estmate over the TESTING set
X_est   = nan(n_bands, n_smp_split, n_drr);
obj_list= cell(1, n_drr);
scores  = nan(1, n_drr);


for k = 1:n_drr    
    % Choose the DRR case for the testing signals
    test_drr  = k;
    aux.vprint(verbose, '--> Testing the *%s* DRR group (%d/%d)...\n', drr.labels{test_drr}, test_drr, n_drr);
    
    % Split for the TESTING data
    X2 = spec_st.Sft{test_drr};
    y2 = squeeze( H(:, test_drr, :) );
    [~, X_test, ~, y_test, splits] = train_test_split(X2, y2, ...
        'n_splits', n_splits, ...    
        'test_grp', test_grp_number ...
    );


    %% RECONSTRUCTION
    % Initialize the reconstruction object
    obj = reconstruct(binwidth, 'f', spec_st.f); 

    
    % Fit the model
    obj.fit(X_train, y_train,...
        'jackknife_flag', jackknife_flag,...
        'n_splits', n_splits, ...
        'fignum', []);
    

    % Predict the spectrogram
    X_est(:,:,k) = obj.predict(y_test);

    scores(k) = goodness(X_test, X_est(:,:,k));
    
    % Save for later analysis
    warning off
    obj_list{k} = struct(obj);
    warning on
    
    
    if verbose
        aux.cprintf('g', '--> RECONSTRUCTION SCORE: %g\n\n', scores(k));
    end  
end


%% Save the data
%{
    'SAVE the analysis!'
    fn.save = sprintf('../.data/reconstruct_SU_bw(%g)_fbands(%d)_win(%g)ms_(%s).mat',...
        binwidth, n_bands, win_size_ms, date);
    
    %obj_st = struct(obj);
    
    save(fn.save, '-v7.3', 'X_test', 'X_est', 'H', 'spec_st', 'obj_list');
    %save(fn.save, 'X_test', 'X_est', 'H', 'spec_st', 'obj');
%}




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


% % Set the color axis to be the same for both spectrograms
% cmin = min([X_test(:); X_est(:)]);
% cmax = max([X_test(:); X_est(:)]);
% caxis(ax(1), [cmin, cmax])
% caxis(ax(2), [cmin, cmax])

linkaxes(ax);


%% Plot #2: count CF
%{
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
%}

% Add ABCD labels
% aux.abc(ax, 48, 'northwestoutside');



%% Plot #3: plot the G filters
obj.plot_reconstruction_filters(spec_st.f, 20+fignum);


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

%% Plot #6: lagged autocorrelation between spectrogram frequency channels & SU 
% norm_cols = 0;
% ac = ac_spec_psth(spec_st.Sft{3}, squeeze(psth_mtx(:,3,:)), norm_cols, 35+fignum, 1e-3*spec_st.f);








