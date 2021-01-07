function [mi, fm, f] = MTF(Sft, varargin)
%
%   function [mi, fm, f] = MTF(Sft, [[f], [binwidth]] or [spec_st])
%
% Input:
%     Sft     : (n_bands x T) spectrogram
%     f       : (n_bands x 1) (Hz) band frequencies of the spectrogram (y-axis)
%     binwidth: (1x1) binwidth of the frequencies
%
%
%
% Description:
%   Calculates the speech transmission index (STI) and the modulation transfer 
%   function (MTF) of the spectrogram (Sft).
%
% Reference:
%     Houtgast, T. and Steeneken, H.J., 1985. A review of the MTF concept in room 
%     acoustics and its use for estimating speech intelligibility in auditoria. 
%     The Journal of the Acoustical Society of America, 77(3), pp.1069-1077.

if isstruct(varargin{1})
    spec_st = varargin{1};
    binwidth= spec_st.binwidth; % (ms)
    f       = spec_st.f;        % (Hz)
    
%     % The spectrogram (Sft) amplitudes are of log scale, so we need to revert it
%     cfloor = spec_st.db_floor;
%     cmax = spec_st.max_Sft;
%     Xft = cmax*10.^((Sft + cfloor)/20);
    
    %Xft     = Sft;
    %warning('--> [STI.m]: cannot calculate the inverse-log amplitudes!!');        
    
    
else    
    binwidth= varargin{1};  % (ms)
    f       = varargin{1};  % (Hz)
    %Xft     = Sft;
    %warning('--> [STI.m]: cannot calculate the inverse-log amplitudes!!'); 
    
end


Xft = Sft;
% Xft = 10.^(Xft/20+spec_st.)*spec_st.max_Sft;

assert(size(Sft,1) == length(f));
[n_bands, nt] = size(Xft);


% Sampling rate of the time-domain axis
fs = 1/(1e-3*binwidth);     % (Hz)

% The 1/3 octave "bandpass filter"
fm   = 2.^((-2:0.2:15)./3);
% fm = logspace(log10(2^-2), log10(2^4), 60);
n_fm = length(fm);

% Modulation transfer function
mi = nan(n_bands, n_fm);

% 2*... : because we are looking at HALF of the modulated frequency domain,
% and because we wish that it would be 1.0 in case of a pure sine wave
mtf = @(x) 2*abs(fft(x)) / sum(x);  

% See: https://github.com/JacobD10/SoundZone_Tools/blob/master/STI.m
freqs = linspace(0, fs/2, nt/2+1);      % (Hz)
freqs(end) = [];    % No nyquist frequency


%% The speech-envelope spectrum
for ii = 1:n_bands
    % Get the ii'th frequency band; because we take the envelope directly 
    % from the spectrogram, we don't need to band-pass (octave wide) filter the
    % input signal and to extract the envelope.
    xj = Xft(ii,:)';
    
    % DEBUG -- MTF normalization 
    % Uncomment these lines to debug the MTF. The normalization of a pure
    % sine wave must be ~1.0 for the sine wave's modulation frequency (fm)  
    % and ~0.0 otherwise.
    %{
        '*** DEBUG ***'
        tt    = linspace(1/fs, size(Sft,2)*1/fs, size(Sft,2))';
        f_idx = 10;
        x_env = 0.5*(1 + cos(2*pi*fm(f_idx)*tt));
        xj    = x_env; 
        %Iavg  = 1/sqrt(2)*mean(xj) * ones(1, n_bands);  
        ff = linspace(1/fs, fs, length(xj))';
        plot(ff, db(abs(fft(xj)))/length(xj));
        xlim([0, fm(end)]);
        xlabel('Frequency (Hz)');
        ylabel('Amp (dB)');
        
        %jj = f_idx;
        %plot(ff, db([abs(fft(xj(:)))/length(xj), abs(fft(filter(B2{jj}, A2{jj}, [1:length(xj)==1])))']))
    %}
    
    % Interpulate instead of the 1/3-octave bandpass
    xf = mtf(xj);
    mi(ii,:) = interp1(freqs, xf(1:floor(end/2)), fm);
end











