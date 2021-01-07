function [y, y_info] = wav(fn, new_sample_rate)
%
% [y, y_info] = load.wav(fn, new_sample_rate)
%
% Input:
%   fn               : (str) filename to load (+ path, if needed)
%   [new_sample_rate]: (Hz,1x1) if given, resample the loaded signal <y> to this
%                       new sampling rate. 
%
% Output:
%   y               : (Nx1) a vector of the signal.%
%   y_info           : (struct) structure that holds information about the
%                     loaded signal, such as, sampling rate, number of 
%                     samples, etc.
%
% Description:
% Loads a WAV file. In case <new_sample_rate> is given, the function will
% resample the WAV file to the desired sampling rate. 
%
%
%

import load.*


if 2 > nargin
    new_sample_rate = [];     % (Hz) the new sample rate
end

assert(~isempty(dir(fn)),...
    sprintf('--> ERROR: can''t find this file:\n\t <%s>!!', fn));


%% Load the big (concatenated) wav file (36 or 40 secs)
[yload, fs] = audioread(fn);
y_info      = audioinfo(fn);

% Down-sample if needed to save space
if ~isempty(new_sample_rate)
    [num, den] = rat(new_sample_rate/fs);
    y = resample(yload, num, den);  % 100k Hz --> 16k Hz
    fs = new_sample_rate;
else
    y = yload;
end

y_info.SampleRate = fs;     % (Hz) update the SampleRate
y_info.TotalSamples = size(y,1);












