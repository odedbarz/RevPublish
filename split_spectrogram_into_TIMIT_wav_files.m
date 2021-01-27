function [split_time_idx, n_splits, tbl_metadata] = ...
    split_spectrogram_into_TIMIT_wav_files(binwidth, duration_sec, margin_sec, Sft)
%
% function [split_time_idx, n_splits, tbl_metadata] = ...
%    split_spectrogram_into_TIMIT_wav_files(duration_sec, [Sft])
%
% Input:
%   binwidth    : (1x1) the time-axis, in milisecond 
%   duration_sec: (1x1) the total duration of the stimulus, in seconds.
%   [margin_sec]: (1x1); default margin_sec, in seconds
%   [Sft]       : (nt x n_bands; optional) a spectrogram to plot; the function uses 
%                 this input for PLOTTING only 
%
% Output:
%   split_time_idx, n_splits, tbl_metadata
%
% Description
% This function uses the TIMIT meta-data table to create indices that
% splits the stimulus into its WAV time intervals.

fs = 1/(1e-3*binwidth);     % (Hz) time axis of the spectrogram

% The total samples in the stimulus
n_time = round(fs*duration_sec);   % (samples)

if 3 > nargin || isempty(margin_sec)
    margin_sec = max(0,floor(0.10 * fs)); % (samples); default margin_sec is 0.5 sec
end

if nargin < 4
    Sft = [];
end

% Load the the metadata file of the stimulus
fn.tbl.path     = load.path_to_data('Stimulus');
fn.tbl.file     = sprintf('metadata_(%d)_wav_(30-Jun-2020)', duration_sec);
dummy           = load( fullfile( fn.tbl.path, fn.tbl.file ) );
tbl_metadata    = dummy.tbl_metadata;
n_splits        = height(tbl_metadata);
% -> convert from samples (fs=16kHz) to samples (fs=1/1e-3*binwidth)
split_time_idx = floor([tbl_metadata.t0, tbl_metadata.t1]/(tbl_metadata.fs(1)/fs));
split_time_idx(1)   = 1;        % start at the beginning of the stimulus
split_time_idx(end) = n_time;   % finish at the end of the stimulus

% Add a small margin from the silence on both side of the chunk
split_time_idx(:,1) = max(1, split_time_idx(:,1)-margin_sec);
split_time_idx(:,2) = min(split_time_idx(:,2)+margin_sec, n_time);


%% DEBUG mode; check that the splitting is OK
if isempty(Sft)
    return;
end

fprintf('-> DEBUG MODE; split_spectrogram_into_TIMIT_wav_files.m:\n');
fprintf('--> split_time_idx:\n');
disp(split_time_idx);


figure(999);
clf;
assert( size(Sft,2) == n_time );
imagesc(Sft); 
colormap winter;
aux.vline(split_time_idx(:,1), 'Color', 'w', 'LineWidth', 3.0);
aux.vline(split_time_idx(:,2), 'Color', 'r', 'LineWidth', 0.5);
for n = 1:n_splits
    t0 = split_time_idx(n,1);
    t1 = split_time_idx(n,2);
    xlim([max(1,t0-40), min(t1+40,n_time)]);
    title(sprintf('Chunk number %d', n));
    pause(1.0);
end

