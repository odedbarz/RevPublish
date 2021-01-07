function h = subsasgn(h, index, val)
% HISTOGRAM/SUBSASGN - Assign values of histogram counts for range of indices
% H.Count(m1:n1,m2:n2) = B assigns the value B to the spike count defined by the indices
%
% For neurograms, other possibilities are H.SyncCount(m:n) = B, 
% H.SpikeCount(m:n) = B, and H.PeriodCount(m:n) = B, 
% and all defined stimulus variables (e.g. H.Frequency(m:n) = B)
%
switch index(1).type
case '.'
    switch index(1).subs
    case 'Count'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h.data(index(2).subs{:}) = val;
       else
           h.data = val;
       end
   case 'SyncCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h.synccount(index(2).subs{:}) = val;
       else
           h.synccount = val;
       end   
   case 'SpikeCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h.spikecount(index(2).subs{:}) = val;
       else
           h.spikecount = val;
       end
   case 'PeriodCount'
       if length(index) >= 2 & strcmp(index(2).type, '()'),
           h.periodcount(index(2).subs{:}) = val;
       else
           h.periodcount = val;
       end
   otherwise
       dim = strmatch(index(1).subs, h.varname, 'exact');
        if length(index) >= 2 & strcmp(index(2).type, '()'),
            h.value{dim}(index(2).subs{:}) = val;
        else
            h.value{dim} = val;
        end
    end
otherwise
    error('Illegal histogram subscripted assignment')
end


