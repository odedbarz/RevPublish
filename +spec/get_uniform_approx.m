function [fstart, fstop, n, err] = get_uniform_approx(f)
%
% This function was taken from MAT:AB's signal processing toolbox, 
%   C:\Program Files\MATLAB\R2017a\toolbox\signal\signal\private\
% 
% I use it specifically for computinf DFT for nonuniform spaced
% frequencies, e.g., log scale frequencies.
%
%
% GETUNIFORMAPPROX get uniform approximation of a frequency vector
%   [Fstart, Fstop, N, Emax] = GETUNIFORMAPPROX(F) returns the first
%   frequency, Fstart, last frequency Fstop, and number of frequency
%   points, N, and the relative maximum error, Emax, between any internal
%   point and a linearly spaced vector with the same number of points over
%   the same range as the input.
%
%   This file is for internal use only and may be removed in a future
%   release.

%   Copyright 1988-2015 The MathWorks, Inc.

fstart = min(f);    % f(1);
fstop  = max(f);    % f(end);
n      = numel(f);
err    = max(abs(f.'-linspace(fstart,fstop,n))./max(abs(f)));
