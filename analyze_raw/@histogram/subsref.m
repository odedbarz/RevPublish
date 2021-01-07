function h = subsref(h, index)
% HISTOGRAM/SUBSREF - Get subsets of a histogram based on indexes
% H(m1:n1,m2:n2) returns a histogram including elements of H defined by
% the indices.  H.Count(m1:n1,m2:n2) returns just the bin counts. 
%
% For neurograms, other possibilities are H.SyncCount(m:n), H.SpikeCount(m:n),
% and H.PeriodCount(m:n), and all defined stimulus variables (e.g. H.Frequency(m:n))
%
switch index(1).type
case '()'
    dmax = ndims(h.data);
    h.data = h.data(index.subs{:});
    
    % adjust offsets
    minsub = [1 1];
    for bwdim = 1:2,
        sub = index.subs{dmax+bwdim-2};
        if isnumeric(sub),
            minsub(bwdim) = min(sub);
        end
    end
    h.offset = h.offset + (minsub-1).*h.binwidth;
    
    % adjust stim values
    if length(h.value) > 0,   % neurogram
        for dim = 1:length(h.value);
            h.value{dim} = h.value{dim}(index.subs{dim});
        end
        h.synccount = h.synccount(index.subs{1:length(h.value)});
        h.spikecount = h.spikecount(index.subs{1:length(h.value)});
        h.periodcount = h.periodcount(index.subs{1:length(h.value)});
    end
    
case '.'
    switch index(1).subs
    case 'Count'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h = h.data(index(2).subs{:});
       else
           h = h.data;
       end
   case 'SyncCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h = h.synccount(index(2).subs{:});
       else
           h = h.synccount;
       end   
   case 'SpikeCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h = h.spikecount(index(2).subs{:});
       else
           h = h.spikecount;
       end
   case 'PeriodCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h = h.periodcount(index(2).subs{:});
       else
           h = h.periodcount;
       end
   otherwise
       dim = strmatch(index(1).subs, h.varname, 'exact');
        if length(index) >= 2 & strcmp(index(2).type, '()'),
            h = h.value{dim}(index(2).subs{:});
        else
            h = h.value{dim};
        end
    end
       
otherwise
    error('Illegal histogram subscripted reference')
end


