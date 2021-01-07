function h = resize(h, new_size)
% RESIZE - Resize histogram or neurogram
% RESIZE(H, new_size) returns a histogram with new size and all counts set to 0
% Bin widths, offsets, spike/sync channels are unchanged
%
dmax = length(new_size);
if dmax < 2,   
    new_size = [1 new_size];    % row histograms are the default
    dmax = 2;
end

if dmax > 2,
    h.data = zeros(new_size);
else    % for efficiency
    h.data = sparse(new_size(1), new_size(2));
end

% resize counts
if dmax > num_dim(h.type),      % neurogram
     csize = new_size(1:dmax-num_dim(h.type));
     if length(csize) < 2, csize = [csize 1]; end
     h.synccount = zeros(csize);
     h.spikecount = zeros([csize(1:end-1) length(h.spikechan)]);
     h.periodcount = zeros(csize);
else
    h.synccount = 0;
    h.spikecount = zeros(size(h.spikechan));
    h.periodcount = 0;
end
 
% resize neurogram variables
for dim = 1:length(h.value),
    if length(h.value{dim}) ~= size(h.data,dim),
        h.value{dim} = [0:size(h.data,dim)-1];
    end
end

