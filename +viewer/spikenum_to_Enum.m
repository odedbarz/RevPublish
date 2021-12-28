function channel_number = spikenum_to_Enum(spikechan, IMPALE_MAX_SUB_CHAN)
%
% function Enum= spikenum_to_Enum(spikechan, IMPALE_MAX_SUB_CHAN)
%
%

if 2 > nargin
    IMPALE_MAX_SUB_CHAN = 6;
end

channel_number = fix(spikechan / IMPALE_MAX_SUB_CHAN) + 1;