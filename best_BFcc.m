function best_BFcc(Sft, y)

Rn = (Sft - mean(Sft,2)) * (y - mean(y));
Rn = (Rn/N)./(std(Sft,[],2)*std(y));   % normalize
[~, kbest] = max(Rn);

% Get the best frequency over all envelopes of the spectrogram
BF(k)= spec_st.f(kbest);    % (Hz)
R(k) = Rn(kbest); %(1,2);             % correlation coefficient of the best frequency  