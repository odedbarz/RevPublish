function h = extract_hist(ng, idx)
% EXTRACT_HIST - Extract histogram from a neurogram
% EXTRACT_HIST(NG, IDX) returns the histogram located at index IDX in neurogram NG
%
if length(ng.value) == 0,
    error('Can''t extract a histogram from another');
end
if length(idx) ~= length(ng.value),
    error('Requested index does not match neurogram variables');
end

% create subscript array
nvar = length(ng.value);
for k = 1:nvar, subs{k} = idx(k); end
for k = 1:num_dim(ng.type), subs{nvar+k} = ':'; end

h = ng;
h.data = shiftdim(ng.data(subs{:}));
if num_dim(h.type) == 1, h.data = h.data.'; end  % restore row vector
h.spikecount = shiftdim(ng.spikecount(subs{1:nvar+1}))';
h.synccount = ng.synccount(subs{1:nvar});
h.periodcount = ng.periodcount(subs{1:nvar});

h.value = {};
h.varname = {};
h.varunit = {};
