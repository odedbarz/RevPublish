function smp = ms2samples(x, Fs)
%
%   function smp = ms2samples(x, Fs)
% 
% Input:
%   x: (1x1) the input is in units of time of [ms]
%
% Output:
%   y: (1x1) the output is in units of time of [seconds]
%
% Description:
%   Converts from [ms] to [sec]
%

smp = floor( aux.ms2sec(x)*Fs );    % floor, not ceil, to truncate extra msecs
smp = max(smp, 1);  % minimum sample in MATLAB is always 1
