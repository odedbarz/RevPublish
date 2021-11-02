function amp_Pa = SPL2Pa(dbspl)
%
%   function amp_Pa = SPL2Pa(amp_dBSPL)
%

amp_Pa = 20e-6 * 10^( dbspl/20 );  


