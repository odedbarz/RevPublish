function coh = interaural_coherence(yrir)
% 
% function coh = interaural_coherence(yrir)
%
% Calculates the interaural coherence. The coherence is the the maximum of
% the absolute value of the cross-correlation of the waveforms reaching the 
% two ears from the source (Lavandier et al., 2008, Hartmann et al., 2005).
% 

len_rir = length(yrir);

coh = nan(len_rir,1);

for kk = 1:len_rir
    y_left_ear = yrir{kk}(:,1);
    y_right_ear = yrir{kk}(:,2);
    
    % The kk'th coherence
    coh(kk) = max( abs(xcorr(y_left_ear, y_right_ear,'coeff')) );
end



