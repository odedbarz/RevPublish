function ng = insert_hist(ng, idx, h)
% INSERT_HIST - Insert histogram into a neurogram
% NG=INSERT_HIST(NG, IDX, H) inserts histogram H at index IDX into neurogram NG.
% Bin values, spike counts and sync counts are updated.
%
if length(h.value) > 0,
    error('Can''t insert a neurogram into another');
end
if length(ng.value) == 0,
    error('Can''t insert a histogram into another');
end
if length(idx) ~= length(ng.value),
    error('Requested index does not match neurogram variables');
end

nd = num_dim(h.type);
nvar = length(ng.value);
ngsize = size(ng.data);
hsize = size(h.data);

if strcmp(ng.type, 'none'),     % initialize neurogram from hist
    
    ng.type = h.type;
    ng.binwidth = h.binwidth;
    ng.offset = h.offset;
    ng.period = h.period;
    ng.order = h.order;
    ng.spikechan = h.spikechan;
    ng.syncchan = h.syncchan;
    if ndims(ng.data) ~= nvar+nd || ...
       any(ngsize(nvar+1:end) ~= hsize(3-nd:2)),
        new_size = [ngsize(1:nvar) hsize(3-nd:2)];
        ng = resize(ng, new_size);
    end
    
else                            % check consistency
    if ~strcmp(h.type, ng.type),
        error('Types do not match');
    end
    if any(ngsize(end-nd+1:end) ~= hsize(end-nd+1:end)),
        error('Sizes do not match');
    end
    % Removed by bard to avoid errors in PH neurograms
%     if any(h.binwidth ~= ng.binwidth),
%         error('Bin widths do not match');
%     end
end


% create subscript array
subs={};
for k = 1:nvar, subs{k} = idx(k); end
for k = 1:nd, subs{nvar+k} = ':'; end

% copy data from histogram to neurogram
ng.data(subs{:}) = full(h.data);
ng.synccount(subs{1:nvar}) = h.synccount;
ng.spikecount(subs{1:nvar+1}) = h.spikecount;
ng.periodcount(subs{1:nvar}) = h.periodcount;


    
    

