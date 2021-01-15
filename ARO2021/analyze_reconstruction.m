%
% analyze_reconstruction.m
%
% Description:
% Compares spectrogram reconstructions. The comparison is between DRY and other
% DRR conditions.
%
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');




%% Load data
% 
%   Name                Size                 Bytes  Class     Attributes
% 
%   H_units          7200x5x10             2880000  double              
%   fn                  1x1                   2206  struct              
%   obj_list            5x12              13504080  cell                
%   sorted_list       241x1                   1928  double              
%   spec_st             1x1               15039029  struct              
%   splits              1x1                  58552  struct              
%   tbl_data          241x20                339094  table               
%
data_type   = 'MUA';       % {'SU', MUA'}
fn_path= '../_data/Reconstruct';
data_type   = upper(data_type);
switch data_type
    case 'SU'
        fn_template = 'reconstruct_SU_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';       
%         unit_list = [10, 25, 50, 103, 150];
        
            unit_list = [10, 25, 50]
        
    case 'MUA'
        fn_template = 'reconstruct_MUA_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';        
        unit_list = [10, 25, 50, 103];  % 150

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end




%% SCORES & STI structures
drr   = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
len_unit_list = length(unit_list);

% Loading first file to get preliminary data
fn_n = sprintf(fn_template, unit_list(1));
aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
warning off
dummy = load( fullfile(fn_path, fn_n), 'splits', 'spec_st', 'stim_st' );
warning on

splits   = dummy.splits;
spec_st  = dummy.spec_st;
stim_st  = dummy.stim_st;
n_splits = splits.n_splits;
n_bands  = spec_st.n_bands;

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)


clear scores  %sti sti2
scores.n_drr        = n_drr;
scores.info         = 'X_dry vs. X_est';
% scores.n_splits         = n_splits;
% scores.n_neuron_list= n_neuron_list;
% scores.n_bands      = n_bands;
scores.n_units      = nan(1, len_unit_list);
scores.CC           = nan(n_drr, n_splits, len_unit_list);
scores.mse          = nan(n_drr, n_splits, len_unit_list);
scores.nmse         = nan(n_drr, n_splits, len_unit_list);
scores.CCf          = nan(n_drr, n_splits, len_unit_list, n_bands);

%{
sti.n_drr        = n_drr;
sti.n_splits         = n_splits;
sti.units        = finfo.n_neurons;
sti.n_neuron_list= n_neuron_list;
sti.test         = nan(n_drr, n_splits);        % X_test is always DRY condition!
sti.test_probe   = nan(n_drr, n_splits);        % X_test is always DRY condition!
sti.est          = nan(n_drr, n_splits, n_neuron_list);
sti.est_probe    = nan(n_drr, n_splits, n_neuron_list);
%}

% For comparing between stimuli of various DRR conditions and the
% reconstructed spectrograms
%{
sti2.n_drr        = n_drr;
sti2.n_splits         = n_splits;
sti2.n_neuron_list= n_neuron_list;
sti2.est          = nan(n_drr, n_splits, n_neuron_list);
%}
scores2           = scores;
scores2.info      = 'X_drr vs. X_est';




%%
for n = 1:len_unit_list
    % Load the file
    fn_n = sprintf(fn_template, unit_list(n));
    aux.cprintf('Comments', '\n-> Loading <%s>...\n', fn_n);
    warning off
    data = load(fullfile(fn_path, fn_n));
    warning on
    
    obj_list= data.obj_list;
    spec_st = data.spec_st;
    %H       = data.H_units;
    splits = data.splits;
    
    % # of neurons
    n_neurons = obj_list{1}.n_neurons;
    assert( unit_list(n) == n_neurons, 'ERROR: something is wrong!' );

    if verbose
        fprintf('--> n_neurons  : %d\n', n_neurons);
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_splits   : %d\n', n_splits);
    end
    

    % Loop over the splits (exclusive intervals)
    scores.n_units(n) = n_neurons;
    for sp = 1:n_splits 
        % Cut out the testing chunk of the spectrogram
        X_dry = spec_st.Sft{drr.dry}(:, sp == splits.idx);    
        
        for k = 1:n_drr    
            rv = drr.ordered(k);
            
            % Cut out the testing chunk of the spectrogram
            X_est                 = obj_list{rv,sp}.X_est;
            gof                   = goodness(X_dry, X_est);                  
            scores.CC(k,sp,n)    = gof.CC;
            scores.mse(k,sp,n)   = gof.mse;
            scores.nmse(k,sp,n)  = gof.nmse;
            scores.CCf(k,sp,n,:) = gof.CCf;            
            %sti.est(rv,sp,n)      = STI(X_est, spec_st);
            %sti.est_probe(rv,sp,n)= STI({X_dry, X_est}, spec_st);
            
            % For comparing between stimuli of various DRR conditions and the
            % reconstructed spectrograms
            X_drr                 = spec_st.Sft{rv}(:, sp == splits.idx);            
            gof2                  = goodness(X_drr, X_est);        
            scores2.CC(k,sp,n)   = gof2.CC;
            scores2.mse(k,sp,n)  = gof2.mse;
            scores2.nmse(k,sp,n) = gof2.nmse;
            scores2.CCf(k,sp,n,:)= gof2.CCf;
            %sti2.est(rv,sp,n)     = STI({X_est, X_drr}, spec_st);
            
            % * # of neurons doesn't make any difference for the STIMULUS
            %if 1 == n
            %    sti.test(rv,sp) = STI(X_drr, spec_st);
            %    sti.test_probe(rv,sp) = STI({X_drr, X_dry}, spec_st);
            %end
            
        end
    end
    

end
 
% Statistics
% Dims: [drr x # splits x # units]
scores.mu.CC    = squeeze( median(scores.CC, 2) );
scores.SE.CC    = squeeze( mad(scores.CC,[],2)/sqrt(n_splits) );
scores.mu.mse   = squeeze( median(scores.mse, 2) );
scores.SE.mse   = squeeze( mad(scores.mse,[],2)/sqrt(n_splits) );
scores.mu.nmse  = squeeze( median(scores.nmse, 2) );
scores.SE.nmse  = squeeze( mad(scores.nmse,[],2)/sqrt(n_splits) );

scores2.mu.CC   = squeeze( median(scores2.CC, 2) );
scores2.SE.CC   = squeeze( mad(scores2.CC,[],2)/sqrt(n_splits) );
scores2.mu.mse  = squeeze( median(scores2.mse, 2) );
scores2.SE.mse  = squeeze( mad(scores2.mse,[],2)/sqrt(n_splits) );
scores2.mu.nmse = squeeze( median(scores2.nmse, 2) );
scores2.SE.nmse = squeeze( mad(scores2.nmse,[],2)/sqrt(n_splits) );




%% BROAD-BAND CC_broadband_envelope
% Calculate the broadband CCs
win_env_lpf_ms = 60;
fs_dwnsmp      = fs_timeaxis;
[R, args]      = broadband_envelopes_CCs(stim_st.Y, stim_st.fs, win_env_lpf_ms, fs_dwnsmp);
CC_broadband_envelope = R(drr.dry, drr.sortby(1:n_drr));





%%















