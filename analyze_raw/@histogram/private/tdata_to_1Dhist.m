function data = tdata_to_1Dhist(h, tdata)
% TDATA_TO_1DHIST - Put time data into 1D histogram
% TDATA_TO_1DHIST(H, TDATA) pust the time data TDATA into the bins of
% the 1D histogram H
% Spike and sync counts are not updated
%
if min(size(h.data)) > 1,
    error('Histogram must be 1-D')
end

if isempty(tdata),
    data = zeros(size(h.data));
    return;
end

dim = which_dim(h);
bins = h.offset(dim) + h.binwidth(dim) * [0:size(h.data,dim)];

n = histc(tdata(:)', bins);
data = n(1:end-1);  % discard last, overflow bin
