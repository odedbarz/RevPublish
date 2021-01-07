function amp_Pa = spl2Pa(dbspl)
%
%   function amp_Pa = spl2Pa(amp_dBSPL)
%

amp_Pa = sqrt(2)* 20e-6 * 10^( dbspl/20 );  % Zilany's code
% amp_Pa = 20e-6 * 10^( dbspl/20 ); 


