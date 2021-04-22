function [Yenv, Y] = filter_bank(x, fs, fc, bw, N)
%
%   [Yenv, Y] = filter_bank(x, fs, [fc], [bw], [N])
%
% Input:
%   x : (Nx1) a time-series to filter
%   fs: (1x1) sampling rate
%   fc: (1xK) center frequencies
%   bw: (1x1 or 1xK) bandwidth of the octave filter-band (i.e., 1/3 or 1)
%   N : (1x1) filter's order 
%
% Output:
%  Yenv: (NxK)
%  Y   : (NxK)
%
% Description:
% The function filters the vector x by octave filters of width bw centered at
% the the frequencies fc.
%
% Default frequencies (Stark & Abeles, 2007)
%     alpha       : 1-13 Hz
%     beta        : 13-30 Hz
%     gamma       : 30-60 Hz
%     higher gamma: 60-100 Hz
%


%% Set the inputs
if 3 > nargin || isempty(fc)
    % Default frequencies (Stark & Abeles, 2007)
    %     alpha       : 1-13 Hz     -> ~6 Hz, bw=1/3
    %     beta        : 13-30 Hz    -> ~16 Hz, bw=1/3
    %     gamma       : 30-60 Hz    -> ~45 Hz, bw=1
    %     higher gamma: 60-100 Hz   -> ~75 Hz, bw=1  
    fc = [6, 16, 45, 75];
end

if 4 > nargin || isempty(bw)
    % bw == 1, octave band; bw == 1/3, 1/3 octave band
    bw = [1/3, 1/3, 1, 1]; 
end

if 5 > nargin || isempty(N)
    N = 3;  % filter's order
end


%%
n_freq = length(fc);
if isscalar(bw)
    bw = bw*ones(1, n_freq);
end

x = x(:);
Y = nan(size(x,1), n_freq);
for k = 1:n_freq       
    [B, A] = octave.design_filter(fc(k), fs, bw(k), N);    
    Y(:,k) = filtfilt(B, A, x);
end
Yenv = envelope(Y);





