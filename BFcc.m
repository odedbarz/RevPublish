function [bf, c, pvalue, kbest] = BFcc(H, Sft, f, debug_mode)
%
% Calculates the best-frequency of the correlation coefficients (BFcc).
%

if nargin < 4
    debug_mode = false
end

N = size(Sft,2);

% extract only the dry measurements 
yi =  squeeze( H );

CC = (Sft - mean(Sft,2)) * (yi - mean(yi));
CC = (CC/N)./(std(Sft,[],2)*std(yi));   % normalize
[~, kbest] = max(CC);

% Get the best frequency over all envelopes of the spectrogram
bf = f(kbest);    % (Hz)
c  = CC(kbest); %(1,2);    

if nargout >= 3
    [~, P] = corrcoef(Sft(kbest,:), yi);
    pvalue = P(1,2);
end

if debug_mode
    plot(f, CC)
    xlabel('Frequency (Hz)');
    ylabel('CC');
    aux.vline(f(kbest));
    title(sprintf('BF: %g', bf));
end