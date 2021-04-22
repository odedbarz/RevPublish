function [B, A] = design_filter(fc, fs, bw, N)
%
%   function [B, A] = design_filter(fc, fs, [bw], [N])
%
% My simple implementation for an octave filter
%
% See: https://stackoverflow.com/questions/35324394/how-to-create-1-3-octave-band-filters-in-matlab-octave
%

if 3 > nargin
    bw = 1;     % set to bw = 1/3 for a 1/3 octave bandpass filter
end

if 4 > nargin || isempty(N)
    N = 3;
end


fl = fc*2^(-bw/2);      % lower cutoffs
fu = fc*2^(+bw/2);      % upper cutoffs

[B, A] = butter(N, [fl, fu]/(fs/2), 'bandpass');