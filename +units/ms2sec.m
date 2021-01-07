function x_sec = ms2sec(x_ms)
%
%   function smp = ms2samples(x_ms)
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

x_sec = 0.001 * x_ms;

