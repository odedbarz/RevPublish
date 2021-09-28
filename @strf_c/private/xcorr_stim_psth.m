function [Ravg, avg, JN] = xcorr_stim_psth(Sft, psth, n_win, varargin)
%
%   function [Ravg, avg, JN] = xcorr_stim_psth(Sft, psth, n_win, varargin)
%
%

%% Set the inputs
p = inputParser;

addRequired(p, 'Sft', @isnumeric);              % (n_bands x n_time) spectrogram
addRequired(p, 'n_win', @isnumeric);            % (1x1, samples) # of lags for xcorr

addOptional(p, 'xcorr_type', 3, @isnumeric);    % (1x1) perform the xcorr in the frequency domain        
addOptional(p, 'jk_flag', 1, @isnumeric);       % (1x1) jack-knight flag
addOptional(p, 'gpu_flag', 0, @isnumeric);      % (1x1) use GPU
addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

parse(p, Sft, n_win, varargin{:});
pars = p.Results;

% xcorr_type = p.Results.xcorr_type;       
% jk_flag    = p.Results.jk_flag;       
% verbose    = p.Results.verbose;

% Don't perform the jackknife if not needed
if 3 > nargout, pars.jk_flag = 0; end



%%
if pars.verbose
    timerVal = tic;
end

% n_chunks = numel( chunk_num );
n_bands  = size(Sft,1);
unbias_single_autocorr_mtx = 0;


%% Remove the means from the frequency bands of the spectrogram
avg.Sft = mean( Sft, 2);
Sft     = Sft - avg.Sft;


%% Remove the PSTH mean
avg.psth = mean(psth);
psth     = psth - avg.psth;


%%

if pars.verbose
    fprintf('--> [autocorr_spectrogram.m]: There are %d segments in this data.\n', n_chunks);
end

% dim1_xcorr = (n_bands+1)*n_bands/2;
Ravg = zeros(n_bands, 2*n_win+1);
avg.bias = zeros(1, 2*n_win+1);


% The autocorrelation matrix
[R, bias] = xcorr_avg( Sft, psth, n_win, pars.xcorr_type, ...
    unbias_single_autocorr_mtx, pars.gpu_flag );

% NOTE:
% The XCORR_AVG computed the <Sft(t),PSTH*(t)>, but we need the
% conjugate transpose of this expression (<Sft*(t),PSTH(t)>). So, using
% Fourier transform properties, we flip the cross-correlation matrix
% along the time domain
R = fliplr( R );

% Aggregate for the average autocorrelation matrix
Ravg = Ravg + R;

% The total correlation bias; this is the factor that the MATLAB's XCORR 
% command + 'unbias' removes
avg.bias = avg.bias + bias;
% end

% Remove bias (along the xcorr window axis)
Ravg = Ravg./avg.bias;

% % Normalize by the # of trials
% Aavg = Aavg/n_trials;

if pars.verbose
    toc(timerVal);
end



%% Jack-knife (JN)
% Calculate CHUNK_NUM cross-correlations. For each k'th JN, calculate the
% cross-correlation without the k'th chunk.
%
if 0 == pars.jk_flag, return; end

if pars.verbose
    fprintf('--> [xcorr_st.m]: staring jack-knife calculation...');
end

% n_chunks= length(chunk_num);
JN.N    = length(chunk_num);
JN.bias = cell(1, n_chunks);
JN.Ravg = cell(1, n_chunks);

for kk = 1:length(chunk_num)
    % The k'th jackknife bias
    JN.bias{kk} = avg.bias - bias{kk};
    
    % The k'th jackknife cross-correlation matrix (stimulus\response 'vector')
    JN.Ravg{kk} = (Ravg.*avg.bias - R{kk}) ./ (eps + JN.bias{kk});
    
end

















