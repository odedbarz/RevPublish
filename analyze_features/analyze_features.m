% test_phonemes.m
%
% Testing the plotting of the phonemes procedures. 




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
splits  = data.splits;
spec_st = data.spec_st;
stim_st = data.stim_st;
% tbl_data= data.tbl_data;
% n_bands = spec_st.n_bands;
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
fn_path_meta = load.path_to_data('Stimulus');
fn_file_meta = sprintf('metadata_(%d)_wav_(30-Jun-2020)', duration_sec);
dummy        = load( fullfile( fn_path_meta, fn_file_meta ) );
tbl_metadata = dummy.tbl_metadata;



%% Loads PYDATA's MAT file (LIBROSA)
% pydata = 
%   struct with fields:
% 
%     stim: [1×1 struct]
%        p: [1×1 struct]
%     harm: [1×1 struct]
%     spec: [1×1 struct]
mat_fn      = 'analyzed_librosa.mat';
mat_full_fn = fullfile( load.path_to_data('data'), 'Analysis', mat_fn);
pydata      = load( mat_full_fn );

% Extract parameters
% pyspec = pydata.spec;
pypitch      = pydata.p;

pyspec      = pydata.spec;
pyspec.f0_smp = 5;
pyspec.f    = pyspec.f(pyspec.f0_smp:end);     % Remove low frequencies
pyspec.Sdb  = pyspec.Sdb(pyspec.f0_smp:end,:);
pyspec.S    = pyspec.S(pyspec.f0_smp:end,:);

% Interpulate to N_TIME
F0          = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pypitch.F0,...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.Sdbi = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.Sdb',...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.Si   = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.S',...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.ti   = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.t,...
    linspace(1,pypitch.duration_seconds,n_time) )';




%% Extract a selected speaker SP
sp = 10;
aux.cprintf('String', '- sp: %d (Speaker number)\n', sp);

% Get the speaker's sex
if contains(tbl_metadata.fn(sp), '_M')
    sp_sex = 'Male';
else
    sp_sex = 'Female';
end

% Replace the SPLITS indices to concur with the TIMIT timings
% idx_sp = sp == splits.idx;      % *REPLACED* indices; time indices for speaker SP
dt_meta = 1/tbl_metadata.fs(sp);
t0 = tbl_metadata.t0(sp) * dt_meta;     % (sec) starting time of the speaker utterance
n0 = 1 + floor(t0/(1e-3*binwidth));	% (smp) starting time of the speaker utterance

t1 = tbl_metadata.t1(sp)/tbl_metadata.fs(sp);   % (sec) ending time of the speaker utterance
n1 = floor(t1/(1e-3*binwidth));                 % (smp) ending time of the speaker utterance

idx_sp = n0:n1;
ti = nonzeros( idx_sp*(1e-3*binwidth) );

% f = spec_st.f;  % (Hz) include all frequencies


% Choose DRR condition:
%   1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}
%drr_k  = 1;     
%
% Spectrograms for the SP speaker
    %Sdry_ = spec_st.Sft{drr.ordered(1)}(:, idx_sp);
    %Sdrr_ = spec_st.Sft{drr.ordered(drr_k)}(:, idx_sp);
    %Sest_ = data.Sest(:, idx_sp);
    %f_ = f;
    
% High resolution Spectrogram (LIBROSA)
Sdry_ = pyspec.Sdbi(:, idx_sp);
f_ = pyspec.f;







%% Checking the PLOT PHONEMES.m function
figh = figure(fignum);
clf;
hz = 1e-3;

% AX #1 -- spectrogram
idx_sub = 4;
subplot(20,1,1+idx_sub:20);

[ax, surf_h] = spec.plot_spectrogram(ti-ti(1), hz*f_, Sdry_,...
    'fontsize', fontsize, 'precision', 2, 'fignum', fignum);

% Add F0s
hold on
plot(ax, ti-ti(1), log2(hz*F0(idx_sp)), 'k', 'LineWidth', 5);
hold off


% AX #2 -- phonemes
ax(2) = subplot(20,1,1:idx_sub);
text_h = plot_phonemes(sp, 'ax', ax(2));

% Add the mean amplitudeto the phoneme's graph
mean_Si = mean(pyspec.Si(:,idx_sp));
hold on
plot(ax(2), ti-ti(1), mean_Si*(0.35/max(mean_Si)),...
    'Color', aux.rpalette('new01'), 'LineWidth', 5);
hold off

drawnow;
pos1 = get(ax(1), 'Position');
ax(2).Position(3) = pos1(3);


% ZOOM-IN
xlim_0 = 0.9;
xlim_1 = 2.2;
xlim([xlim_0, xlim_1]);





%%
figh = figure(10+fignum);
clf;
hz = 1e-3;

















