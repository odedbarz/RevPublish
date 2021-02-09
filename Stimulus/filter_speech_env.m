function [yenv, yoct, hoct] = filter_speech_env(yt, varargin)
%
%   function yenv = filter_stim_env(yt, varargin)
%
% Description:
% Calculates the envelope of yt using an octave-band filter.
%

%% Parse the input
p = inputParser;

addRequired(p, 'yt', @(x) isvector(x) && isnumeric(x));

addOptional(p, 'N', 5, @isscalar);
addOptional(p, 'foct', 1500);        % (Hz)
addOptional(p, 'Fs', 100e3);        % (Hz)
addOptional(p, 'fignum', 0);        % plot the resulting envelope filter

parse(p, yt, varargin{:});

N      = p.Results.N;            
foct   = p.Results.foct;    
Fs     = p.Results.Fs;            
fignum = p.Results.fignum;            


%%
[B, A] = octave.octdsgn(foct, Fs, N);
yoct = filtfilt(B, A, yt);
yenv = envelope( yoct );

% Lowpass Butterworth filter
Hd = filt.LPF_BW;
yenv = Hd.filter(yenv);

% % Remove the LPF's lag
% [~, lpf_lag_idx] = max(Hd.filter(1:1000==1));
% yenv = circshift(yenv, -lpf_lag_idx);

hoct = [];
if 3 <= nargout 
    % Get the impulse response, if desired
    hoct = filter(B, A, [1==1:length(yt)]);
end

%% Plot (debug mode)
if isempty(fignum) || (0 == fignum), return; end
figure(99); clf;
subplot(2,1,1);
len_yt = length(yt);
t = 1e3*linspace(0, 1/Fs*len_yt, len_yt)';
plot(t, [yt, yoct, yenv]);
xlabel('Time (msec)');
ylabel('Amp.');
legend('$y(t)$', '$y_{oct}(t)$', '$y_{env}(t)$');
%
subplot(2,1,2);
f = linspace(0, Fs, len_yt)';
plot(f, 20*log10(abs(fft([yt, yoct, yenv]))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');
xlim([0, max(1000,min(Fs/2,2*foct))]);













