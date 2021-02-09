function s = NoiseBB_stim(mod_freq,band_cent,pass_band,varargin)
% Create a binaural beating noise, based on the algorithm in Akeroyd
% (2010).
% Inputs:
%   - mod_freq = Modulation frequency (Hz)
%   - band_cent = Center frequency of the unshifted passband (Hz)
%   - pass_band = Bandwidth of the passband for the noise

%%% Note (10/19/14) -- Akeroyd's method uses bandpass noise, currently I'm using
%%% broadband noise

% Variables
Fs = 100000; % sampling rate
t_end = 1; % duration
tok = randi(100); % randomization token

% Parse varargin
if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

%% Create stimulus
rng(tok); % randomize generator

% Create band-pass noise
pass_inds = pass_band*t_end; % number of indexes in the pass band
cent_ind = band_cent*t_end;
RE = randn(pass_inds+1,1);
IM = randn(pass_inds+1,1);

% Apply band-pass noise in one ear, place frequency-shifted version in the other
S = zeros(t_end*Fs,2);
L_inds = cent_ind-round(pass_inds/2):cent_ind+round(pass_inds/2); % original
R_inds = cent_ind-round(pass_inds/2)+mod_freq*t_end:cent_ind+round(pass_inds/2)+mod_freq*t_end; % frequency shifted
S(L_inds,1) = RE + 1i*IM;
S(R_inds,2) = RE + 1i*IM;
S = [S(1:t_end*Fs/2,:); flipud(S(1:t_end*Fs/2,:))];
s = real(ifft(S));

% Normalize so max(abs(s)) = 1
s = s./(max(max(abs(s))));

%%% Using a method similar to Siveke et al (2008)
% % Create Gaussian noise
% bb = randn(t_end*Fs,1);
% L_chan = bb;
% 
% % Frequency shift the components in the noise
% BB = fft(bb);
% cBB = circshift(BB(1:floor(length(BB)/2)),modfreq*t_end);
% R_chan = real(ifft([cBB; flipud(cBB)]));
% 
% % Adjust the R_chan to have the same RMS as L_chan
% rmsL = sqrt(mean(L_chan.^2));
% rmsR = sqrt(mean(R_chan.^2));
% R_chan = R_chan*(rmsL/rmsR);

% Save sound
% s = [L_chan R_chan];