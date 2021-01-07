%
% test_ASD.m
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

 
 
 %% Define the file information
finfo.type       = 'MUA';           	% {'PSTH', 'MUA'}
finfo.n_neurons  = [10 25 50 100 150];  	% # of units to load from existing files
finfo.binwidth   = 10;               	% (ms)
finfo.n_bands    = 30;
finfo.win_size_ms= 5;
finfo.n_splits   = 12;
finfo.date       = '09-Jun-2020';        % '18-May-2020'
finfo.path       = '../.data/reconstruct/';
fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(%d)ms_fbands(%d)_win(%d)ms_splits(%d).mat';
n_neuron_list    = length(finfo.n_neurons);


% === Load first batch to get preliminary info
finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(1), finfo.binwidth, ...
    finfo.n_bands, finfo.win_size_ms, finfo.n_splits);
dummy = load(fullfile(finfo.path, finfo.fn), 'obj_list');
% ===

% D1: # of DRR cases checked 
% D2: # of  simulations 
[n_drr, n_JK] = size(dummy.obj_list);
n_splits      = n_JK;
n_bands       = dummy.obj_list{1}.n_bands;

% clear scores dummy
% scores.mu = nan(n_drr, n_neuron_list);
% scores.SE = nan(n_drr, n_neuron_list);
scores.n_units = nan(1, n_neuron_list);
scores.CC      = nan(n_drr, n_JK, n_neuron_list);
scores.mse     = nan(n_drr, n_JK, n_neuron_list);
scores.nmse    = nan(n_drr, n_JK, n_neuron_list);
scores.CCf     = nan(n_drr, n_JK, n_neuron_list, n_bands);


%%
for n = 1 %:n_neuron_list
    % Load the data
    %          H: [3600×6×150 double]
    %      X_est: [30×180×5×20 double]
    %     X_test: [30×180 double]
    %   obj_list: {5×20 cell}
    %    spec_st: [1×1 struct]
    % tbl_impale: [437×19 table]
    finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(n), finfo.binwidth, ...
        finfo.n_bands, finfo.win_size_ms, finfo.n_splits);

    aux.cprintf('g', '\n-> Loading <%s>...\n', finfo.fn);
    data = load(fullfile(finfo.path, finfo.fn));
    obj_list= data.obj_list;
    spec_st = data.spec_st;
    H       = data.H;
    splits  = data.splits;
    
    % # of neurons
    n_neurons = obj_list{1}.n_neurons;

    % Make sure that the loaded data and the requested data have the same
    % parameters
    binwidth = spec_st.binwidth;     % (ms)
    assert(binwidth == finfo.binwidth);

    n_bands = spec_st.n_bands;
    assert(n_bands == finfo.n_bands);

    win_size_ms = spec_st.win_size_ms; 
    assert(win_size_ms == finfo.win_size_ms);

    if verbose
        fprintf('\n');
        fprintf('--> n_neurons  : %d\n', finfo.n_neurons(n));
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_JK       : %d\n', n_JK);
    end
    

    %% Get the test groups
    assert( n_splits == size(H,1)/size(data.X_est,2) );
    assert( n_JK == splits.n_grps );

    


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

    
    
    

    %%
    scores.n_units(n) = n_neurons;
    for k = 1:n_drr
        for m = 1:n_JK
            % Cut out the testing chunk of the spectrogram
            % * n_JK == test_grp_number 
            X_test = spec_st.Sft{obj_list{1}.train_drr}(:, m == splits.idx);    
            X_est  = obj_list{k,m}.X_est;
            gof    = goodness(X_test, X_est);
        
            scores.CC(k,m,n)    = gof.CC;
            scores.mse(k,m,n)   = gof.mse;
            scores.nmse(k,m,n)  = gof.nmse;
            scores.CCf(k,m,n,:) = gof.CCf;
            
            %{
            'DEBUG'
            imagesc([X_test; X_est]);
            
            imagesc([zca(X_test')'; zca(X_est')']);
            
            plot([X_test(15,:)', X_est(15,:)']);
    
            plot(1e-3*spec_st.f, squeeze(scores.CCf(drr.sortby(1:5),1,1,:)), '.-')
            
            a = @(n,m) squeeze(scores.CCf(drr.sortby(1:5),n,m,:));
            a_ = @(n,m) squeeze(scores.CCf(drr.sortby(1),n,m,:));
            b = @(n,m) a(n,m)'./a_(n,m);
            plot(1e-3*spec_st.f, b(3,1), '.-');
            %}
        end
    end
    

end







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















