function [strf, strf_noncausal] = calc_strf(obj, X_train, r_train, varargin)
%
%   function strf = calc_strf(obj, X_train, r_train, [gamma, xcorr_type, verbose, fignum])
%
% Description:
% Calculates one STRF given a training set (X_train) and a response batch (r_train).
%


%% 
p = inputParser;

% === Required inputs ===
addRequired(p, 'X_train', @isnumeric);              % (n_bands x n_time) spectrogram power
addRequired(p, 'r_train', @isnumeric);         % (n time x 1) psth of the response

% Tolerance for the ridge regression
addOptional(p, 'gamma', obj.ridge.gamma, @isnumeric);      % (Mx1)       

% Algorithm type to calculate the STRF 
addOptional(p, 'algo_type', obj.algo_type, @isstr);      % % {'ASD', 'regression'}             

% If 1 (true) then output a causal STRF(s), that is, all STRFs are zero for
%   t<0 ms. Otherwise, if 0, return the "non-causal" STRF(s).
addOptional(p, 'iscausal', obj.iscausal, @isnumeric);      % (1x1)       

% Choose the cross-correlation type to use
addOptional(p, 'xcorr_type', obj.xcorr_type, @isnumeric);    % (1x1) correlation methos       

% For the plottings ONLY
addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

parse(p, X_train, r_train, varargin{:});
pars            = p.Results;


% Update the object, if needed
obj.iscausal    = pars.iscausal;
obj.ridge.gamma = pars.gamma;
obj.xcorr_type  = pars.xcorr_type;
obj.algo_type   = pars.algo_type;

% For the plots\debug
verbose = pars.verbose;
fignum  = pars.fignum;



obj.X_train = X_train;
obj.r_train = r_train;



%%
if strcmpi('regression', obj.algo_type)
    obj.asdstats = [];

    % N_LAGS is always odd
    assert( 0 == rem(obj.n_lags+1, 2), '--> [calc_strf.m]: N_LAGS must always be ODD!');

    % Set the lags according to the causality flag
    if obj.iscausal
        n_win_half = obj.n_lags-1;
        n_win = 2*(obj.n_lags-1) + 1;
    else
        n_win_half = (obj.n_lags-1)/2;
        n_win = obj.n_lags;        
    end
    assert( n_win_half == fix(n_win_half), '--> [STF_C->CALC_STRF.m]: n_win_half must be an integer!!' );

    obj.Cavg = autocorr_stim(X_train,...
        n_win_half,...      the size of the temporal window (# of lags)
        'xcorr_type', obj.xcorr_type);

    % Smoothing + FFT to the AC matrix & CCs
    assert(size(obj.Cavg, 2) == n_win);

    % Smoothing window for the FFT
    win_smooth = hanning( n_win );  

    % Spectrogram's autocorrelation (of frequency bands)
    perform_fftshift = 1;        
    obj.Cw = corr_to_fft(obj.Cavg, win_smooth, perform_fftshift);

    % Spectrogram & stimulus cross-correlation
    [obj.Ravg, obj.bias] = xcorr_stim_psth(X_train, r_train,...
        n_win_half,...      the size of the temporal window (# of lags)
        'xcorr_type', obj.xcorr_type );

    % FFT of the correlation averaged matrix
    perform_fftshift = 0;    
    obj.Rw = corr_to_fft(obj.Ravg, win_smooth, perform_fftshift);

    % Inverse the correlation matrix Ridge regression 
    [obj.strf, obj.strf_noncausal] = decorrelation(obj.Cw, obj.Rw, obj.ridge.gamma, obj.iscausal, verbose);

    
elseif strcmpi('ASD', obj.algo_type)
    obj.ridge = [];
    
    Xcnv    = convmtx_spec(X_train, obj.lags);
    nks     = [obj.n_bands, obj.n_lags];  % number of filter pixels along [cols, rows]
    minlens = [2; 2];  % minimum length scale along each dimension

    [kasd, obj.asdstats] = fastASD(Xcnv', r_train, nks, minlens);
    obj.strf = reshape(kasd, obj.n_bands, []);

    
else
    error('--> [strf_c/calc_strf.m]: unrecognized ALGO_TYPE (%s) parameter!!!', pars.algo_type);
    
end



%% Set the output:
if 1 <= nargout
    strf = obj.strf;   
    strf_noncausal = obj.strf_noncausal;
end


%% DEBUG Mode
if ~isempty(fignum)
    figure(fignum);
    clf;
    add_margins = 1;
    obj.plot_strf(add_margins);
end







