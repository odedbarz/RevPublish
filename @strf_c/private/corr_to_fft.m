function Cw = corr_to_fft(Ct, win_smooth, perform_fftshift)
%
%   function Xw = corr_to_fft(Xt, win_smooth, perform_fftshift)
%

if 3 > nargin || isempty(perform_fftshift)
    perform_fftshift = 1;
end
 

n_lags = size(Ct, 2);

% Temporal smoothing window before applying FFT
if ~exist('win_smooth','var') || isempty(win_smooth)
    win_smooth = hanning(n_lags);
end
win_smooth = win_smooth(:)';    % win_smooth must be a row vector

% Smooths (LPF) the boundaries so they will gradually reduce to zero. 
Xsm = Ct .* win_smooth;

% The XCORR (time domain) produces a 'centered' reponse (peak at the center
% of the vcector); thus, when converting to the frequency domain, this 
% introduces a linear phase. To remove this linear phase, we remove the 
% temporal shift in the time domain (i.e., using IFFTSHIFT).
if perform_fftshift
    Xsh = ifftshift( Xsm, 2 );
else
    Xsh = Xsm;
end

% DFT for each row (temporal domain\lags)
Cw = fft(Xsh,[], 2);


