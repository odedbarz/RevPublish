function [b,a] = SynchronyFilter(Fco, Fs)
% SYNCHRONYFILTER - Create filter matching AN synchony data
% Uses a 4th-order all-pole lowpass filter
% usage: [b,a] = SynchronyFilter(Fco, Fs)
% Fco    3-dB Filter cutoff frequency (default 1000 Hz)
% Fs     Sampling rate in Hz (default 20000)
% b, a   Filter coefficients
%
if nargin < 1, Fco = 1000; end
if nargin < 2, Fs = 20000; end

n = 4;    % filter order
Fco = Fco / sqrt(2^(1/n)-1);  % correct for cumulation of cutoff attenuation
zp = exp(-2*pi*Fco/Fs);     % pole
a0 = [1, -zp];   % first-order lowpass

b = (1-zp)^n;     % give unit gain at DC
a = a0;
for k = 1:n-1, a = conv(a, a0); end
