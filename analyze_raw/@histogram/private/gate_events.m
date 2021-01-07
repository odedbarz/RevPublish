function [t, ch] = gate_events(t, ch, gate, syncchan)
% GATE_EVENTS - Gate time events
% Usage: [tg, chg] = gate_events(t, ch, gate, syncchan)
%   t, ch       raw events times and channels
%   gate        structure containing 3 fields: 'type,', 'width', and 'delay'
%               'type' is either 'none', 'time' (default), 'PST', or 'trial'
%               'time' selects events over a range of recording times
%               'PST' selects over a range of peristimulus times
%               'trial' selects a range of stimulus presentations
%               'delay', 'width'  start and duration of gate in msec (except for 'trial')
%   syncchan    sync pulse channel for 'PST' and 'trial' gates
%   tg, chg     gated events
%
switch gate.type
case 'none'
    return;
case 'time'
    k = find(t >= gate.delay & t < gate.delay+gate.width);
case 'PST'
    if nargin < 4, syncchan = 0; end
    tps = pstimes(t, ch, syncchan);
    k = find(ch == syncchan | (tps >= gate.delay & tps < gate.delay+gate.width));
case 'trial'
    if nargin < 4, syncchan = 0; end
    trialno = cumsum(ch == syncchan);
    k = find(trialno >= gate.delay & trialno < gate.delay+gate.width);
otherwise
    error(['Unsupported gate type: ' type]);
end
    
t = t(k);
ch = ch(k);
    