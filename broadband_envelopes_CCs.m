function [R, args] = broadband_envelopes_CCs(Y, fs, win_env_lpf_ms, fs_dwn)
%
%   [R, args] = function broadband_envelopes_CCs(Y, fs, [win_env_lpf_ms], [fs_dwn])
%
% Input:
% Y             : (nt x n_drr): stimuli for various DRR conditions
% fs            : (1x1; Hz) sample frequency of the stimuli Y 
% win_env_lpf_ms: (1x1) size of the low-pass window (hanning) for th
%                 envelope
% fs_dwn        : (1x1; Hz) frequency downsample to extract the envelope
%
%
% Output:
% R   : (n_drr x n_drr) correlation coefficients between stimuli in Y
% args: (struct) all the arguments that was used to extract the broadband
%       envelopes
%
%
% Description: 
% Broadband envelopes coefficient correlations (CCs).
%

if 3 > nargin
    % downsample window
    win_env_lpf_ms = 60;   % (ms)
end

if 4 > nargin    
    fs_dwn = 1000;  % (Hz) default LP frequency
end


%%
% Window
args.win_env_lpf_ms = win_env_lpf_ms;                     % (ms)
args.n_win = ceil( 1e-3*args.win_env_lpf_ms * fs);     % # of samples for the window
args.win_type = 'hanning';
args.lpf_win = eval(sprintf('%s(%d)', args.win_type, args.n_win));

% Rectify the matrix with all the stimuli
Y_rectify = max(0, Y);

% Create the envelope; filter the stimuli with the window
Yenv = filtfilt(args.lpf_win, 1, Y_rectify);

% Resampling ratio
[num, den] = rat(fs_dwn/fs);
args.Yenv_dwn = resample(Yenv, num, den);

[R, args.P] = corrcoef(args.Yenv_dwn);


