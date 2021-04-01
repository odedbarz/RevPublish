%
% analyze_setup.m
%
% Description:
% Loads analyzed data to the workspace. These data are used by various scripts
% for further analysis and plotting.
%
%


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
data_type   = 'SU';       % {'SU', MUA'}
switch data_type
    case 'SU'
        n_units = 103;  %[10, 25, 50, 103]; 
        
    case 'MUA'
        n_units = 103;   %[10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end


fn_path = load.path_to_data('Analysis');
fn_name = sprintf('analyzed_cc_%s_units(%d).mat', data_type, n_units);
fn_fullfile = fullfile( fn_path, fn_name );
warning off
data = load(fn_fullfile);
warning on

aux.cprintf('String', '\nLoading analyzed data...\n');
aux.cprintf('String', 'data_type: %s\n', data_type);
aux.cprintf('String', 'n_units  : %d\n', n_units);
aux.cprintf('String', 'filename : %s\n', fn_name);
aux.cprintf('String', 'path     : %s\n\n', fn_path);

fprintf('Loading the CC arrays...\n');


% Print this to get more information about the data
if exist('verbose','var') && verbose
    aux.cprintf('Comments', 'Data Info:\n');
    aux.cprintf('Comments', data.info);
end

% Extract data: 
CC      = data.tbl.CC.Variables;         % Sdry vs. Sest, averaged over time
CC2     = data.tbl.CC2.Variables;        % Sdrr vs. Sest
CC3     = data.tbl.CC3.Variables;        % Sdry vs. Sdrr

CCt     = data.CCt;         % compares dry vs. est, as a function of time
CCt2    = data.CCt2;        % compares drr vs. est
CCt3    = data.CCt3;        % compares drr vs. est

% non-nans indices
idx_nn = ~isnan(CCt(:,5)) & ~isnan( CCt3(:,5) );


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
voiced_probs_i = pypitch.voiced_probs_i;
voiced_flag_i  = pypitch.voiced_flag_i;

% Get all voiced intervals
idx_F0i = ~isnan(F0i) .* F0i>0;

% make it a boolean vector
idx_vc = 1  == voiced_flag_i;  

% Convenient for plotting
vc_nans = double(idx_vc); 
vc_nans(1 ~= vc_nans) = nan;

%{
sp     = 1;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
idx_sp = idx_fun(sp);      % indices; time indices for speaker SP
drr_k  = 5;      % 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}
%}

 
%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( spec_st );







