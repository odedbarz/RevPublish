function [y] = rect_upsample(x,upsamp)
% Upsamples the 1-D signal x by zero-padding the fourier transform of the
% signal to the sampling rate which is 'upsamp' times the original sampling
% rate.  If x has multiple columns, the upsampling and interpolation is
% applied to each columns separately
% Inputs:
%   - x = signal
%   - upsamp = the ratio of the desired sampling rate to the original
% Outputs:
%   - s = the upsampled signal

cols = size(x,2);
X = fft(x);
Y = [X(1:round(length(X)/2),cols); zeros(length(X)*(upsamp-1),cols); X(round(length(X)/2)+1:end,cols)];
y = real(ifft(Y))*upsamp;