% test_phonemes.m
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
% CC      = data.CC;         % compares dry vs. est
% CC2     = data.CC2;        % compares drr vs. est
% CCt     = data.CCt;         % compares dry vs. est
% CCt2    = data.CCt2;        % compares drr vs. est
splits  = data.splits;
spec_st = data.spec_st;
stim_st = data.stim_st;
% tbl_data= data.tbl_data;
% n_bands = spec_st.n_bands;
% n_time  = length(splits.idx);


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






%% Extract a selected speaker SP
sp     = 10;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP

% Choose DRR condition:
%   1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}
drr_k  = 1;      

% Get the speaker's sex
if contains(tbl_metadata.fn(sp), '_M')
    sp_sex = 'Male';
else
    sp_sex = 'Female';
end

% Replace the SPLITS to concur with the TIMIT timings
dt_meta = 1/tbl_metadata.fs(sp);
t0 = tbl_metadata.t0(sp) * dt_meta;     % (sec) starting time of the speaker utterance
n0 = 1+1 + floor(t0/(1e-3*binwidth));	% (smp) starting time of the speaker utterance

t1 = tbl_metadata.t1(sp)/tbl_metadata.fs(sp);   % (sec) ending time of the speaker utterance
n1 = floor(t1/(1e-3*binwidth));                 % (smp) ending time of the speaker utterance

idx_sp = n0:n1;
ti = idx_sp*(1e-3*binwidth);

f = spec_st.f;  % (Hz) include all frequencies

% Spectrograms for the SP speaker
Sdry_ = spec_st.Sft{drr.ordered(1)}(:, idx_sp);
Sdrr_ = spec_st.Sft{drr.ordered(drr_k)}(:, idx_sp);
Sest_ = data.Sest(:, idx_sp);




% Plot Spectrogram + phonemes
figh = figure(fignum);
clf;
hz = 1e-3;

% AX #1
idx_sub = 4;
ax = subplot(20,1,1+idx_sub:20);

[ax, surf_h] = spec.plot_spectrogram(ti-ti(1), hz*f, Sdry_,...
    'fontsize', fontsize, 'precision', 2, 'fignum', fignum);
% ax = gca;
% img = imagesc(t_, hz*f_, Sdry_ );
% set(ax, 'YDir', 'normal');
% set(ax, 'FontSize', fontsize);
% arrayfun(@(X) colormap(X, 'jet'), ax );


% New table of phonemes of a selected speaker
tbl_phonemes = [table(tbl_metadata.phn{sp}, 'VariableNames', {'phn'}),...
    table(tbl_metadata.phn_start{sp} - tbl_metadata.t0_timit(sp), 'VariableNames', {'phn_start'}),...
    table(tbl_metadata.phn_end{sp} - tbl_metadata.t0_timit(sp), 'VariableNames', {'phn_end'}),...
];

tbl_phonemes = [tbl_phonemes, ...
    table(tbl_phonemes.phn_start * 1/tbl_metadata.fs(sp), 'VariableNames', {'t0'}),...
    table(tbl_phonemes.phn_end * 1/tbl_metadata.fs(sp), 'VariableNames', {'t1'})...
];

duration_sp_sec = (tbl_metadata.t1(sp) - tbl_metadata.t0(sp))/tbl_metadata.fs(sp);

% AX #2
ax(2) = subplot(20,1,1:idx_sub);
set(ax(2), 'XTickLabels', '');
set(ax(2), 'YTickLabels', '');

aux.vline(tbl_phonemes.t0(1), 'ax', ax(2), 'Color', 'k');
aux.vline(tbl_phonemes.t1, 'ax', ax(2), 'Color', 'k');
linkaxes(ax, 'x');

d_times = tbl_phonemes.t1 - tbl_phonemes.t0;
t_phn = tbl_phonemes.t0 + 0.5*d_times;
phn_for_text = aux.mName2latex(tbl_phonemes.phn);
text(ax(2), t_phn, mean(ylim)*ones(length(t_phn),1), phn_for_text,...
    'HorizontalAlignment', 'center', 'FontSize', fontsize);

% xlim([0.5, 2.5]);

drawnow;
pos1 = get(ax(1), 'Position');
ax(2).Position(3) = pos1(3);































