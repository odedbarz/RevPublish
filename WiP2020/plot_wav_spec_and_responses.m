% plot_wav_spec_and_responses.m
%
% Plots stimulus, spectrogram, and responses for the presentation.

clc
fignum = 11;
verbose = 1;
addpath('../');
FigSetup;



%%
fn.load.path    = '../.data';
fn.load.file    = 'data_MUA_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
fn.load.fullfile= fullfile( fn.load.path, fn.load.file );
data            = load(fn.load.fullfile);
spec_st         = data.spec_st;
duration_sec    = 1e-3*spec_st.duration_ms;  % (sec) stimulus duration to use
        
% Load the stimuli
dummy   = load(fullfile('../.data/stimulus', 'data_stimuli_duration(36_and_40)sec.mat'));
stim_st = dummy.stim_list{1};
assert(spec_st.duration_ms == stim_st.duration_ms);

drr = get_DRR_list_and_indices;
%%
sp = 6;     % speaker # to plot

[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(spec_st.binwidth, duration_sec);

% What is the speaker's sex?
if contains(tbl_metadata.fn{sp}, '_M')
    test_sex = 'Male';
else
    test_sex = 'Female';
end

% Spectrogram
spec_times     = split_time_idx(sp,:);
spec_idx       = spec_times(1):spec_times(end);
t_spec         = spec_st.t(spec_idx) - spec_st.t(spec_idx(1));

% Stimulus
stimulus_times = spec_idx .* (stim_st.fs/(1/(spec_st.binwidth*1e-3)));
stim_idx       = stimulus_times(1):stimulus_times(end);
t_stim         = stim_st.t(stim_idx) - stim_st.t(stim_idx(1));


%% Plot stimulus WAVE
add_envelope = true;
fontsize = 34;
drr_case = 1;
drr_idx = drr.ordered(drr_case);

figure(2*(drr_case-1) + fignum);
clf;

y = stim_st.Y(stim_idx, drr_idx);
plot(t_stim, y);
title( sprintf('"%s"\n(%s Speaker)', tbl_metadata.txt{sp}(9:end-2), test_sex) );
set(gca, 'YTickLabel', '');

% Add ENVELOPE
if add_envelope
    yenv = envelope(y, 200, 'rms');
    hold on
    plot(t_stim, 2.0*yenv, 'LineWidth', 6);
    hold off
end
xlabel('Time (sec)');    

display_name_1 = sprintf('$s(t)$ (%s)', drr.labels{drr_idx});
display_name_2 = sprintf('envelope (rms)');
legend(display_name_1, display_name_2);
set(gca, 'FontSize', fontsize);

% set FIGURE position
set(gcf, 'Position', [66, 247, 1852, 731]);



%% Plot stimulus SPECTROGRAM
figure(10+fignum);
clf;
Xft = spec_st.Sft{drr.dry}(:,spec_idx);

    Xft = spec_st.Sft{2}(:,spec_idx);
    Xft = max(0, Xft-30);

nolabels = 0;
spec.plot_spectrogram(t_spec, 1e-3*spec_st.f, Xft, 10+fignum, 0, fontsize); 
ylabel('Frequency (kHz)');
title(sprintf('Speaker %d', sp));






