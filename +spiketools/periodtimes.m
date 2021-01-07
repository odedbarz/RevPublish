function [tp, nperiod] = periodtimes(t, ch, period, syncchan)
% PERIODTIMES - Get event times modulo stimulus period
% PERIODTIMES(T< CH, PERIOD) returns the event times modulo PERIOD in msec
% with an origin at the beginning of recording (T=0).
% PERIODTIMES(T< CH, PERIOD, SYNCCHAN) returns the event times modulo PERIOD 
% with an origin reset at each event in SYNCCHAN.
%
if nargin < 3,
    error('Specifiy event times, channels and period');
end
if nargin < 4 | isempty(syncchan),
    tp = rem(t, period);
    nperiod = (t(end)-t(1))/period;
else
    tp = rem(pstimes(t, ch, syncchan), period);
    ksync = find(ch == syncchan);
    if length(ksync) > 1,
        nperiod = length(ksync)*mean(diff(t(ksync)))/period;
    else
        nperiod = (t(end)-t(ksync))/period;
    end
end
