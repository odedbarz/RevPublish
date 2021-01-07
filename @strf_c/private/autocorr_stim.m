function [Cavg, avg_Sft, xcorr_bias] = autocorr_stim(Sft, n_win, varargin)
%
%   [Cavg, avg_Sft, xcorr_bias] = autocorr_stim(Sft, n_win, varargin)
%
%

%% Set the inputs
p = inputParser;

addRequired(p, 'Sft', @isnumeric);              % (n_bands x n_time) spectrogram
addRequired(p, 'n_win', @isnumeric);            % (1x1, samples) # of lags for xcorr

addOptional(p, 'xcorr_type', 3, @isnumeric);    % (1x1) perform the xcorr in the frequency domain     
addOptional(p, 'gpu_flag', 0, @isnumeric);      % (1x1) use GPU
addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

parse(p, Sft, n_win, varargin{:});
pars = p.Results;

% xcorr_type  = p.Results.xcorr_type;       
% verbose     = p.Results.verbose;



%%
if pars.verbose
    timerVal = tic;
end

n_bands = size(Sft,1);
unbias_single_autocorr_mtx = 0;

% Calculate the overall mean along the frequency bands
avg_Sft = mean( Sft, 2);

% Remove the frequency means from the spectrogram
Sft = Sft - avg_Sft;

Cavg = zeros((n_bands+1)*n_bands/2, 2*n_win+1);
xcorr_bias = zeros(1, 2*n_win+1);
   
% The k'th autocorrelation matrix
[Ck, bias_k] = xcorr_avg( Sft, [], n_win, pars.xcorr_type, ...
    unbias_single_autocorr_mtx, pars.gpu_flag );

% Aggregate for the average autocorrelation matrix
Cavg = Cavg + Ck;

% The correlation bias
xcorr_bias = xcorr_bias + bias_k;

% Remove bias (along the xcorr window axis)
Cavg = Cavg./xcorr_bias;

% % Normalize by the # of trials
% Aavg = Aavg/n_trials;

if pars.verbose
    toc(timerVal);
end









