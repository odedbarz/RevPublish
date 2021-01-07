function amp_dBSPL = Pa2SPL( amp_Pa )
%
%   function amp_dBSPL = Pa2SPL( amp_Pa )
%
% According to:
%   amp_Pa = 20e-6 * 10^( amp_dBSPL/20 ); 
%

amp_dBSPL = 20*log10(amp_Pa/20e-6);