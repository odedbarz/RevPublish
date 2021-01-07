function count = count_events(ch,chan)
% count = count_events(ch,chan)
% counts events on specified channel(s)
% ch can be either a cell array or numeric
% if ch is a cell array, chan can be either a scalar or a vector with the
% same length as ch
% if ch is numeric, chan can be a vector of any size

if iscell(ch)
    
    %check for valid input
    if length(chan) == 1
        chans = chan*ones(length(ch),1);
    elseif length(chan) == length(ch)
        chans = chan;
    else
        error('For cell array input, chan must be either a scalar or a vector of the same length as ch!')
    end
    
    count = zeros(length(ch),1);
    %count events
    for k=1:length(ch)
        count(k) = sum(ch{k} == chans(k));
    end
else
    count = zeros(length(chan),1);
    for k = 1:length(chan),
        count(k) = sum(ch == chan(k));    
    end
end