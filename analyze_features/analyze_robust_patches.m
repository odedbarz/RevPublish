%
% analyze_robust_patches.m
%
% Description:
% Compares spectrogram reconstructions. The comparison is between DRY and other
% DRR conditions.
%
%

clc
fignum = 10;
verbose = 1;

setup_environment('../');




%% Plot properties
fontsize = 32;
fontsize_big = 42;
fontsize_bigger = 64;

markersize = 24;




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
        %unit_list = [10, 25, 50, 103];
        unit_list = 103;
        
    case 'MUA'
        fn_template = 'reconstruct_MUA_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';        
        %unit_list = [10, 25, 50, 103, 150, 240];
        unit_list = 241

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end




%% Initialization
drr = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
drr_idx = drr.sortby(1:n_drr);
% len_unit_list = length(unit_list);


% Loading first file to get preliminary metadata
%        H_sorted: [7200×5×10 double]
%              fn: [1×1 struct]
%        obj_list: {5×12 cell}
%     sorted_list: [1×241 double]
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [241×20 table]
fn_n = sprintf(fn_template, unit_list(1));
aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
warning off
dummy = load( fullfile(fn_path, fn_n), 'splits', 'spec_st', 'stim_st', 'tbl_data' );
warning on

splits   = dummy.splits;
spec_st  = dummy.spec_st;
stim_st  = dummy.stim_st;
tbl_data = dummy.tbl_data;
n_splits = splits.n_splits;
n_bands  = spec_st.n_bands;
n_units  = height(tbl_data);

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
t           = spec_st.t;            % (sec)


%%
lag_list_smp = 1    %[1, 5, 10];      % samples
len_lag_list = length(lag_list_smp);

clear scores scores2
scores.info    = 'dry vs. est';
scores.n_drr   = n_drr;
scores.n_units = [];
scores.CC      = nan(n_drr, n_splits, len_lag_list);
scores.CC      = nan(n_drr, n_splits, len_lag_list);

scores2 = scores;
scores2.info = 'dry vs. est';


for n = 1:len_lag_list
    n_lags = lag_list_smp(n);
    
    % Load the file
    fn_n = sprintf(fn_template, unit_list(n));
    aux.cprintf('Comments', '\n-> Loading <%s>...\n', fn_n);
    if ~exist('data', 'var')
        warning off
        data = load(fullfile(fn_path, fn_n));
        warning on
        obj_list = data.obj_list;
        spec_st  = data.spec_st;
        splits   = data.splits;
    end
    
    % # of neurons
    assert( unit_list(n) ==  obj_list{1}.n_neurons, 'ERROR: something is wrong!' );
    n_neurons = obj_list{1}.n_neurons;

    if verbose
        fprintf('--> n_neurons  : %d\n', n_neurons);
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_splits   : %d\n', n_splits);
        fprintf('--> lags_ms    : %g ms\n', n_lags * binwidth);
    end
    

    % Loop over the splits (exclusive intervals)
    scores.n_units(n) = n_neurons;
    for sp = 1:n_splits 
        % Cut out the testing chunk of the spectrogram
        Sdry = spec_st.Sft{drr.dry}(:, sp == splits.idx);    
        Adry  = zca(im2col(Sdry, [n_bands, n_lags], 'sliding'));

        for k = 1:n_drr    
            rv = drr.ordered(k);
            
            % Cut out the testing chunk of the spectrogram
            Sdrr = spec_st.Sft{rv}(:, sp == splits.idx);            
            Adrr = zca(im2col(Sdrr, [n_bands, n_lags], 'sliding'));     % apply ZCA for CC

            Sest = obj_list{rv,sp}.X_est;
            Aest = zca(im2col(Sest, [n_bands, n_lags], 'sliding'));     % apply ZCA for CC
            
            % CCs:
            cc_all_dry_est     = corrcoef(Sdry(:), Sest(:));  
            scores.CC(k,sp,n)  = cc_all_dry_est(1,2);
            
            cc_all_drr_est     = corrcoef(Sdrr(:), Sest(:));   
            scores2.CC(k,sp,n) = cc_all_drr_est(1,2);
            
            plot(t, [mean(Adry .* Aest)', mean(Adrr.*Aest)'] )
            

        end
    end
    
    
    

end
 













