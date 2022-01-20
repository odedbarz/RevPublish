function [bf, r] = BFcc(H, Sft, f)
%
% Calculates the best-frequency of the correlation coefficients (BFcc).
%

N = size(Sft,2);

% extract only the dry measurements 
yi =  squeeze( H );

Rn = (Sft - mean(Sft,2)) * (yi - mean(yi));
Rn = (Rn/N)./(std(Sft,[],2)*std(yi));   % normalize
[~, kbest] = max(Rn);

% Get the best frequency over all envelopes of the spectrogram
bf = f(kbest);    % (Hz)
r  = Rn(kbest); %(1,2);    