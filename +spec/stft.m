function [Sx, f, t, spec_st] = stft(x, fs, varargin)
%
%   function [Sx, f, t, spec] = stft(x, fs, ...)
%
% Input:
%   x           : (Nx1) the stimulus
%   fs          : (1x1) sample rate of the signal x
%   [n_bands]   : (1x1) number of frequency samples; also the window's size
%   [binwidth]  : (1x1) bin-width along the time axis
%   [fignum]    : (1x1, int) if empty plot nothing; if int, create a new 
%                 figure with this number
%
%
% Output:
%   t & f: time (sec) and frequency (Hz) axes of the spectrogram.
%   Sf   : (PxQ) ONE-SIDED power spectrogram (without the phase)
%   spec  : (struct) a structure that contains the two-sided spectrogram and
%          all parameters that have been used for creating the spectrogram.
%
% Description:
%   Computes the spectrogram of x.
%


%% Check the Input
p = inputParser;

addRequired(p, 'x', @isvector);
addRequired(p, 'fs', @isnumeric);

% # frequencies
addOptional(p, 'f', [], @(x) isnumeric(x));

% # of frequency bands (and size of the window)
addOptional(p, 'n_bands', 61, @(x) isnumeric(x));

% The temporal window 
addOptional(p, 'win', [], @isnumeric);

% # of frequency bands (and size of the window); if WIN is defined, it
% overwrites this parameter
addOptional(p, 'n_win', [], @isnumeric);

% (logical) remove the mean from each frequency band (that is, mean of time axis for each band)
addOptional(p, 'remove_mean', 1, @(x) isnumeric(x) || islogical(x));

% (logical) retains only the power (abs); removes the phase
addOptional(p, 'abs_spectrum', 1, @(x) isnumeric(x) || islogical(x));   

% (Hz) the spectrogram's time axis (x-axis) final sample rate
addOptional(p, 'binwidth', 1, @isnumeric);      % BINWIDTH instead of WIN OVERLAP

% (Hz) highest frequency in the spectrogram
addOptional(p, 'highfreq', 8000, @isnumeric);   

% (Hz) lowest frequency in the spectrogram
addOptional(p, 'lowfreq', 250, @isnumeric);    

% lower dB for threshold    
addOptional(p, 'thr_db', 80, @isnumeric);                  

% (logical) output the LOG of the spectrogram
addOptional(p, 'do_log', 1, @(x) isnumeric(x) || islogical(x));                  

% Plot 
addOptional(p, 'fignum', 0, @isnumeric);                  

parse(p, x, fs, varargin{:});


f           = p.Results.f;            
n_bands     = p.Results.n_bands;            
remove_mean = p.Results.remove_mean;            
abs_spectrum= p.Results.abs_spectrum;          
binwidth    = p.Results.binwidth;     
highfreq    = p.Results.highfreq;       
lowfreq     = p.Results.lowfreq;       
thr_db      = p.Results.thr_db;       
do_log      = p.Results.do_log;       
win         = p.Results.win;   
n_win       = p.Results.n_win;   
fignum      = p.Results.fignum;  

% x MUST be a column vector
if isrow(x), x = x(:); end   


%% Check & set the input variables
% frequency step size; this parameter sets the step size along the time
% domain
fs_step = 1/units.ms2sec(binwidth);     

% # of time steps at binwidth Hz
assert(floor(fs/fs_step) == fs/fs_step,...
    '--> Choose a new fs_time or add a code to fix this issue (by resampling x)!');
        
if ~isempty(f)
    % Ignoring the N_BANDS input!    
    n_bands = length(f);
else
    f = linspace(lowfreq, highfreq, n_bands)';
end

% Define the window
if isempty(win)
    % Ignoring the N_WIN input!    
    win = hanning(n_win);
else
    n_win = length(win);
end

% The step size, as a function of the BINWIDTH, along the time-axis
increment = floor(fs/fs_step);  
n_overlap = n_win - increment;



%% Start calculating the spectrogram 
% Initiate the spectrogram's matrix
% tiled_mtx = buffer(x, win_size, win_size-increment);
[tiled_mtx, t] = spec.get_STFT_columns(x, n_win, n_overlap, fs);

% Apply the temporal window to each column
tiled_win_mtx = win .* tiled_mtx;


%% DFT; two-sided
dim = 1;
Sx = spec.compute_DFT(tiled_win_mtx, f, fs, dim);

% Keep the spectrum amplitude, without the phase
if true == abs_spectrum
    Sx = abs(Sx);     
end

% Get the LOG of spectrogram
if ~isempty(thr_db) && true == do_log
    spec_power = max(abs(Sx(:)));  % maximum value (log10(1)==0)
    Sx = max(0, 20*log10(Sx/spec_power) + thr_db );  % (Sx/spec_power) is already in dB

elseif true == do_log
    Sx = 20*log10(Sx);

end    

% Remove the mean from each frequency band?
if true == remove_mean
    Sx = Sx - mean(Sx,2);
end



%% OUTPUT; set the output parameters
if 3 < nargout
    spec_st.Sx = Sx;          % the two-sided spectrogram
    spec_st.t = t;              % (sec) time axis
    spec_st.f = f;              % (Hz) frequency axis
    
    spec_st.n_bands        = size(Sx,1);              
    spec_st.remove_mean    = remove_mean;             
    spec_st.half_freq      = half_freq;
    spec_st.abs_spectrum   = abs_spectrum;      
    spec_st.binwidth       = binwidth;   
    spec_st.highfreq       = highfreq;      
    spec_st.lowfreq        = lowfreq; 
    spec_st.thr_db         = thr_db;          
    spec_st.duration_ms    = 1e3*spec_st.t(end);
    spec_st.fs             = 1/diff(spec_st.t(1:2));      % (Hz)    
end



if isempty(fignum) || 0 == fignum, return; end     % return? or plot the result?
%% Plot
figure(fignum);
% clf;

imagesc(t, 1e-3*f, Sx);
set(gca, 'YDir', 'normal');
colormap jet
xlabel('Time (sec)');
ylabel('Freaquency (kHz)');
if 1 == do_log
    title('One-sided Log Power Spectrogram');
else
    title('One-sided Power Spectrogram');
end


