function h = clear_hist(h)
% CLEAR_HIST - Clear histogram
% CLEAR_HIST(H) returns a histogram in which all data and counts are set to 0
% Size, bin widths, offsets, channels, stimulus variables, etc.. are unchanged
%
dsize = size(h.data);
if ndims(h.data) == 2, 
    h.data = sparse(size(h.data,1),size(h.data,2));
else
    h.data = zeros(size(h.data));
end
h.synccount = zeros(size(h.synccount));
h.spikecount = zeros(size(h.spikechan));
h.periodcount = zeros(size(h.periodcount));

