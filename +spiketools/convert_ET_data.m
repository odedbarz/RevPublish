function data = convert_ET_data(t, ch, order)
% CONVERT_ET_DATA(T, CH, ORDER) converts sun-style ET data to PC-style ET
% data for sequence of stimuli.  T and CH are the times and channel arrays.
% ORDER is the order vector (nonnegative values only)
%
if any(order < 0), 
    error('Negative order values - select completed measurements'); 
end
order = order + 1; % stored order vector starts at 0, not 1

change = find(ch == 2); % Channel 2 signals stimulus change
if length(change) < 2*max(order), 
    error('Incompatible order and channel data');
end
change = [change; length(t)+1];  % to avoid indexing problems

kstart = change(2*order) + 1;
kend = change(2*order+1) - 1;
for k = 1:length(order),
    krange = kstart(k):kend(k);
    data.t{k} = t(krange);
    data.ch{k} = ch(krange);
end