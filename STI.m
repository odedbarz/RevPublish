function [sti, mi, fm] = STI(Sft, spec_st)
%
%   function [sti, midx, fm] = STI(Sft, spec_st)
%
% Input:
%     Sft     : (n_bands x T) spectrogram
%     f       : (n_bands x 1) (Hz) band frequencies of the spectrogram (y-axis)
%     binwidth: (1x1) binwidth of the frequencies
%
%
%
% Description:
%   Calculates the speech transmission index (STI) and the modulation transfer 
%   function (MTF) of the spectrogram (Sft).
%
% Reference:
%     Houtgast, T. and Steeneken, H.J., 1985. A review of the MTF concept in room 
%     acoustics and its use for estimating speech intelligibility in auditoria. 
%     The Journal of the Acoustical Society of America, 77(3), pp.1069-1077.



if iscell(Sft)
    % From dB to amp.
    Sft{1} = 10.^((Sft{1}+spec_st.db_floor)/20)*spec_st.max_Sft;
    %Sft{1} = Sft{1}(spec_st.f<=4e3,:);    % frequencies: 250Hz - 4kHz
    assert( 0 == nnz(isnan(Sft{1}(:))) );
    
    Sft{2} = 10.^((Sft{2}+spec_st.db_floor)/20)*spec_st.max_Sft;
    %Sft{2} = Sft{2}(spec_st.f<=4e3,:);    % frequencies: 250Hz - 4kHz
    %spec_st.f = spec_st.f(spec_st.f<=4e3);
    assert( 0 == nnz(isnan(Sft{2}(:))) );
    
    % Use a PROBE as a referance
    [mi, fm, f] = MTF(Sft{1}, spec_st);
    mi_probe = MTF(Sft{2}, spec_st);
    mi = mi./(eps + mi_probe);
    
else
    % From dB to amp.
    Sft = 10.^((Sft+spec_st.db_floor)/20)*spec_st.max_Sft;
    %Sft = Sft(spec_st.f<=4e3,:);    % frequencies: 250Hz - 4kHz
    %spec_st.f = spec_st.f(spec_st.f<=4e3);
    assert( 0 == nnz(isnan(Sft(:))) );
    
    % No Probe
    [mi, fm, f] = MTF(Sft, spec_st);
    
end




%% STI; Speech transmission index
% Following Houtgast et al., 1985; 

% 1. Transformation to apparent S/N ratio
% SNR = pow2db(mi);
% 
% or...
% Theoretically (Houtgast et al., 1985), 0.0 <= midx <= 1.0
mi_ = min(1, max(0.0 + eps, mi));    	% [!!] practically it doesn't!
SNR = 10*log10(mi_./(1 - mi_));         % (decibels) Eq. 2 in Houtgast et al., 1985
assert( any(~isnan(sum(SNR(:)))) );

% 2. Limiting (S/N)' to a 30-dB range
SNR = min(15, max(-15, SNR));
assert(all(all(SNR<=15 | SNR>=-15)));

% 3. Octave-band-specific mean (S/N)'
band_audibility = mean(SNR, 2);    % (n_bands x 1)

% 4. overall mean <(S/N)'>
% Interpulates using the weights from  Steeneken et al., 1980.
%
% See Steeneken et al., 1980, A physical method for measuring speech-
%   transmission quality
sti_freq    = [0.125, 0.250,   0.5,     1,     2,     4,     8];    % (kHz)
sti_weights = [ 0.13,  0.14,   0.11, 0.12,  0.19,   0.17, 0.14];    % weight (paper & STI standard)
wf          = interp1(sti_freq, sti_weights, 1e-3*f);
wf          = wf/sum(wf);
assert(all(~isnan(wf)));

% re-weighting Houtgas et al. 1985 weights (Goldsworthy & Greenberg, 2004)
% %{
alpha = 1.0;             % no re-weighting
%n_fm = length(fm);      % # of modulated frequencies (1/3 octave bandpass)
%n_f  = length(wf);      % # of bands used for the MTF
%alpha = n_f*linspace(1, n_fm, n_f)'/n_fm;

% Normalize the weights so that the STI will be on the [0,1] interval
wf = (alpha .* wf)/sum(alpha .* wf);
%}

SNavg = wf' * band_audibility;

% 5. Conversion to STI
sti = ( SNavg + 15 )/30;









