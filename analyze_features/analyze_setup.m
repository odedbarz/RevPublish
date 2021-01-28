%
% analyze_cc.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions.
%
%

clc
fignum = 10;
verbose = 1;

setup_environment('../');



%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 5;



%% Load the data
%              CC: [7200×5 double] 	% compares dry vs. est
%             CC2: [7200×5 double]	% compares drr vs. est
%             CCt: [7200×5 double]  % compares dry vs. est along the time domain
%            CCt2: [7200×5 double]  % compares drr vs. est along the time domain
%             PPt: [7200×5 double]
%            PPt2: [7200×5 double]
%     patch_width: 1
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [103×20 table]
data_type   = 'MUA';       % {'SU', MUA'}
switch data_type
    case 'SU'
        n_units = 103;
        
    case 'MUA'
        n_units = 241;   %[10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end


% patch_width: # of samples along the x-axis (time) for each patch to use for comparison
patch_width = 1;    


fprintf('Loading the CC arrays...\n');
fn_path = '../_data/Analysis/';
fn_name = sprintf('analyzed_cc_%s_units(%d)_patchWidth(%d).mat', data_type, n_units, patch_width);
fn_fullfile = fullfile( fn_path, fn_name );
warning off
data = load(fn_fullfile);
warning on

% Print this to get more information about the data
if verbose
    aux.cprintf('Comments', 'Data Info:\n');
    aux.cprintf('Comments', data.info);
end

% Extract data: 
CC      = data.CC;         % compares dry vs. est
CC2     = data.CC2;        % compares drr vs. est
CCt     = data.CCt;         % compares dry vs. est
CCt2    = data.CCt2;        % compares drr vs. est
% splits  = data.splits;
spec_st = data.spec_st;
stim_st = data.stim_st;
tbl_data= data.tbl_data;
n_bands = spec_st.n_bands;
n_time  = spec_st.n_time;


% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
t           = spec_st.t;            % (sec)
drr         = get_DRR_list_and_indices; 
n_drr       = drr.n_drr;
drr_labels  = drr.labels(drr.ordered);



%% Load the the META-DATA of the stimuli (TIMIT) files
duration_sec = 1e-3*data.stim_st.duration_ms;
% fn_path_meta = load.path_to_data('Stimulus');
% fn_file_meta = sprintf('metadata_(%d)_wav_(30-Jun-2020)', duration_sec);
% dummy        = load( fullfile( fn_path_meta, fn_file_meta ) );
% tbl_metadata = dummy.tbl_metadata;
[split_times_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, duration_sec, 0);

idx_fun = @(SP) split_times_idx(SP,1):split_times_idx(SP,2);



%% Loads PYDATA's MAT file (LIBROSA)
[pyspec, pypitch, pydata] = load_pydata(n_time);
F0i = pypitch.F0i;

% Get all voiced intervals
idx_F0i = ~isnan(F0i) .* F0i>0;




 










