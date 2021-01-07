function [Sx, f, t, spec_st] = multitaper(x, fs, varargin)
%
%   function [Sx, f, t, spec] = multitaper(x, fs, n_fft, params, fignum)
%
% Input:
%   x       : (Nx1) the stimulus
%   fs      : (1x1) sample rate of the signal x
%   n_bands : (1x1) number of frequency samples; also the window's size
%   binwidth: (1x1) bin-width along the time axis
%   [params]: (struct)  
%     * binwidth = 1;    % (ms) the spectrogram's time axis (x-axis) final sample rate
%     * highfreq = 8e3;  % (Hz) highest frequency in the spectrogram
%     * lowfreq = 200;   % (Hz) lowest frequency in the spectrogram
%     * thr_db = 80;     % (dB) lower threshold relative to (normalized to 1) peak of the spectrogram
%     * do_log = 1;      % (logical) output the LOG of the spectrogram
%   [fignum]:                   % (1x1, int) if empty plot nothing; if int, create a new figure with this number

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
% Reference:
%   https://www.youtube.com/watch?v=qIlVpZB4u8A


%% Check the Input
p = inputParser;

addRequired(p, 'x', @isvector);
addRequired(p, 'fs', @isnumeric);

% # frequencies
addOptional(p, 'f', [], @(x) isnumeric(x));

% # of frequency bands (and size of the window)
addOptional(p, 'n_win', 61, @(x) isnumeric(x));


addOptional(p, 'nw', 2, @(x) isnumeric(x));

% (logical) remove the mean from each frequency band (that is, mean of time axis for each band)
addOptional(p, 'remove_mean', 1, @(x) isnumeric(x) || islogical(x));

% (Hz) the spectrogram's time axis (x-axis) final sample rate
addOptional(p, 'binwidth', 1, @isnumeric);      % BINWIDTH instead of WIN OVERLAP

% lower dB for threshold    
addOptional(p, 'thr_db', 80, @isnumeric);                  

% (logical) output the LOG of the spectrogram
addOptional(p, 'do_log', 1, @(x) isnumeric(x) || islogical(x));                  

% Plot 
addOptional(p, 'fignum', 0, @isnumeric);                  

parse(p, x, fs, varargin{:});


f           = p.Results.f;            
nw          = p.Results.nw; 
remove_mean = p.Results.remove_mean;            
n_win       = p.Results.n_win;            
binwidth    = p.Results.binwidth;     
thr_db      = p.Results.thr_db;       
do_log      = p.Results.do_log;       
fignum      = p.Results.fignum;  

% x MUST be a column vector
if isrow(x), x = x(:); end   


%% Check & set the input variables
fs_time = 1/units.ms2sec(binwidth);

% # of time steps at binwidth Hz
assert(floor(fs/fs_time) == fs/fs_time,...
    '--> Choose a new fs_time or add a code to fix this issue (by resampling x)!');

% The step size, as a function of the BINWIDTH, along the time-axis
increment = floor(fs/fs_time);  
n_overlap = n_win - increment;
% n_steps   = floor( (length(x) - increment)/n_overlap );


%% Start calculating the spectrogram 
% Initiate the spectrogram's matrix
% tiled_mtx = buffer(x, win_size, win_size-increment);
[tiled_mtx, t] = spec.get_STFT_columns(x, n_win, n_overlap, fs);


%% DFT; two-sided
%dim = 1;
%Sx = spec.compute_DFT(tiled_win_mtx, f, fs, dim);

[Sx, f_] = pmtm(tiled_mtx, nw, f, fs);
if ~isempty(f) && norm(f-f_) > 5*eps
    warning('--> [spec.multitaper]: output frequencies are not as desired!!');
else 
    f = f_;
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
    
    spec_st.nw             = nw;              
    spec_st.remove_mean    = remove_mean;             
    %spec_st.half_freq      = half_freq;
    %spec_st.abs_spectrum   = abs_spectrum;      
    spec_st.binwidth       = binwidth;   
    %spec_st.highfreq       = highfreq;      
    %spec_st.lowfreq        = lowfreq; 
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


