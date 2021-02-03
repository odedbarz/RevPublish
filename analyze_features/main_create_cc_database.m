%
% main_create_cc_database.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions.
%
%

close all;
clc;
clear all;


fignum = 10;
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
data_type   = 'SU';       % {'SU', MUA'}
fn_path   = '../_data/Reconstruct';
data_type = upper(data_type);
switch data_type
    case 'SU'
        fn_template = 'reconstruct_SU_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';       
        n_units = 50;  	% [10, 25, 50, 103]
        
    case 'MUA'
        fn_template = 'reconstruct_MUA_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';        
        n_units = 241;   %[10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end




%% Initialization
drr = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
drr_idx = drr.sortby(1:n_drr);


% Loading first file to get preliminary metadata
%        H_sorted: [7200×5×10 double]
%              fn: [1×1 struct]
%        obj_list: {5×12 cell}
%     sorted_list: [1×241 double]
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [241×20 table]
fn_n = sprintf(fn_template, n_units);
aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
warning off
dummy = load( fullfile(fn_path, fn_n), 'splits', 'spec_st', 'stim_st', 'tbl_data' );
warning on

% Extract data: 
splits   = dummy.splits;
spec_st  = dummy.spec_st;
stim_st  = dummy.stim_st;
tbl_data = dummy.tbl_data;
n_splits = splits.n_splits;
n_bands  = spec_st.n_bands;
n_time   = length(splits.idx);
% n_lags: # of samples along the x-axis (time) for each patch to use for comparison

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
t           = spec_st.t;            % (sec)


%%
% Pearson correlation for average of ALL the split
CC = nan(n_drr, n_splits);          % Sdry-vs-Sest
CC2 = nan(n_drr, n_splits);         % Sdrr-vs-Sest
CC3 = nan(n_drr, n_splits);         % Sdry-vs-Sdrr

% Instant Pearson correlation as a function o DRR
CCt = nan(n_time, n_drr);    	% Sdry vs Sest
PPt = nan(n_time, n_drr);

CCt2 = nan(n_time, n_drr);     	% Sdrr vs Sest
PPt2 = nan(n_time, n_drr);

CCt3 = nan(n_time, n_drr);   	% Sdry vs Sdrr ;only STIMULI
PPt3 = nan(n_time, n_drr);



% Concatenate all the estimations into one big matrix 
Sest = nan(n_bands, n_time, n_drr);
assert( 0 == sum(sum(isnan(spec_st.Sft{1}))), '--> There are NaNs in the DRY spectrograms!' );


% Load the file
fn_n = sprintf(fn_template, n_units);
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
assert( n_units ==  obj_list{1}.n_neurons, 'ERROR: something is wrong!' );
n_neurons = obj_list{1}.n_neurons;

if verbose
    fprintf('\n--> data_type  : %s\n', data_type);
    fprintf('--> n_neurons  : %d\n', n_neurons);
    fprintf('--> binwidth   : %d\n', binwidth);
    fprintf('--> n_bands    : %d\n', n_bands);
    fprintf('--> win_size_ms: %d\n', win_size_ms);
    fprintf('--> n_drr      : %d\n', n_drr);
    fprintf('--> n_splits   : %d\n', n_splits);
end


debug_flag = 0;

% Loop over the splits (exclusive intervals)
for sp = 1:n_splits 
    idx_sp = sp == splits.idx;

    % Cut out the testing chunk of the spectrogram
    Sdry_ = spec_st.Sft{drr.dry}(:, idx_sp);    

    for k = 1:n_drr    
        rv = drr.ordered(k);

        % Cut out the testing chunk of the spectrogram
        Sdrr_ = spec_st.Sft{rv}(:, idx_sp);           

        Sest_ = obj_list{rv,sp}.X_est;

        % Concatenate all reconstructions\estimations into one big 3D
        Sest(:, idx_sp, k) = Sest_;
        
        % AVERAGED Pearson correlation coefficient:
        % DRY-vs-EST
        gof = goodness(Sdry_(:), Sest_(:));  
        CC(k,sp) = gof.CC;
        CCf(k,sp) = gof.CC;
        
        % DRR-vs-EST
        gof2 = goodness(Sdrr_(:), Sest_(:));   
        CC2(k,sp) = gof2.CC;
        
        % DRR-vs-EST
        gof2 = goodness(Sdry_(:), Sdrr_(:));   
        CC3(k,sp) = gof2.CC;

        % INSTANTANEOUS Pearson correlation coefficient:
        %CCt(idx_sp,k) = mean(Adry .* Aest)';
        [CCt(idx_sp,k),  PPt(idx_sp,k)]  = corrcoef_array(Sdry_, Sest_);
        [CCt2(idx_sp,k), PPt2(idx_sp,k)] = corrcoef_array(Sdrr_, Sest_);
        [CCt3(idx_sp,k), PPt3(idx_sp,k)] = corrcoef_array(Sdry_, Sdrr_);

        % DEBUG:
        if debug_flag
            figure(99);
            clf;
            tidx = t(idx_sp);
            plot(tidx, [CCt(idx_sp,k), CCt2(idx_sp,k)] );
            legend('CC(Dry-est)', sprintf('CC(%s-est)', drr.labels{rv}));
            xlabel('Time (sec)');
            ylabel('CC');
            title(sprintf('%s %d', data_type, n_units));
        end
    end
end


assert( 0 == sum(isnan(Sest(:))), ...
    '--> There are NaNs in one (or more) of the reconstructed spectrograms!' );

info = ['- All DRR dimensions are already SORTED!\n',...
    '- Sest: (n_bands, n_time, n_drr)\n',...
    '\n',...
    '- CC : Sdry-vs-Sest, averaged over time\n',...
    '- CC2: Sdrr-vs-Sest, averaged over time\n',...
    '- CC3: Sdry-vs-Sdrr, averaged over time\n',...
    '- CCt : Sdry-vs-Sest, as a function of time\n',...
    '- CCt2: Sdrr-vs-Sest, as a function of time\n',...
    '- CCt3: Sdry-vs-Sdrr, as a function of time\n'...
];




    
%% Save the analysis results
% %{
fprintf('\n- About to save the analysis results...\n');
fn_path = '../_data/Analysis/';
fn_name = sprintf('analyzed_cc_%s_units(%d).mat', data_type, n_units);
fn_fullfile = fullfile( fn_path, fn_name );

% Save the results for that 
save(fn_fullfile, ... 
    ... '-v7.3', ...
    'info',...              general info about the data
    'splits', ...           saves the chunks\intervals\speakers
    'stim_st', ...          stimulus data
    'spec_st', ...          spectrogram's structue with all relevant data
    'tbl_data', ...         a table with all neurons in the data set            
    'Sest',...              reconstructed\estimated spectrograms for all DRRs             
    'CC', 'CCt', 'PPt',...
    'CC2', 'CCt2', 'PPt2', ...
    'CC3', 'CCt3', 'PPt3' ...
); 

fprintf('--> File <%s> SAVED!\n', fn_fullfile);
%}
 













