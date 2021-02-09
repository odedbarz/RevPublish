function [Yrir, yrir] = ConvRIR(Ydry, rir, norm_max)
%
%   function yrir = ConvRIR(y, hrir)
%
% Description:
%   Convolve the room\reverberant impulse response (RIR) with the signal Y 
%

if 3 > nargin
    norm_max = true;
end

% Pad with zeros before the filtering
if iscell(Ydry)
    assert( prod(size(Ydry) == size(rir)) );
    Yrir = cellfun(@(H,YD) fftfilt([H; zeros(size(H))], YD), Ydry, rir, 'UniformOutput', false);
else
    Yrir = cellfun(@(H) fftfilt([H; zeros(size(H))], Ydry), rir, 'UniformOutput', false);
end

yrir = [];
if (1 == size(Yrir,1)) && (1 == size(Yrir,1))
    yrir = Yrir{1};   
end


%% Normalization; avoid clipping when saving to a WAV file
if norm_max
    sig_norm = @(x) x./max(abs(x));    
    Yrir = cellfun(@(YY) sig_norm(YY), Yrir, 'UniformOutput', false);
end