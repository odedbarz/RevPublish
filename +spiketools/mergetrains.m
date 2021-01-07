function [t, ch] = mergetrains(varargin)
% MERGETRAINS - Merge spike trains from different neurons to form labeled event train
% [T, CH] = mergetrain(T1, T2, T3, ...) takes the individual spike trains
% T1, T2, T3, and merge them into a sorted event list, where T are the event
% times, and CH the channels
%
t = varargin{1};
ch = ones(size(t));

for k = 2:length(varargin),
    t = [t; varargin{k}];
    ch = [ch; repmat(k, length(varargin{k}), 1)];
end

[t, ks] = sort(t);
ch = ch(ks);

    