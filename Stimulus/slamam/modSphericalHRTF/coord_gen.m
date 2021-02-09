function coord = coord_gen(azimuth,elev,dist,reference)
%
% COORD_GEN convert azimuth, elevation, distance coordinates to 
%	    Cartesian coordinates.
%
%    COORD = coord_gen(AZIMUTH,ELEV,DIST,REFERENCE) converts the set of
%    M coordinates specified in the M by 1 vectors AZIMUTH, ELEV, and DIST
%    into an M by 3 matrix, COORD, of the [X Y Z] cartesian coordinates.  
%    The coordinates are assumed to be stated relative to the 1 by 3
%    Cartesian reference point, REFERENCE (default [0 0 0]).  
%    DIST and REFERENCE are assumed to be in the same units.  
%    AZIMUTH and ELEV are in degrees.
%
%    For this coordinate system, the XY plane is the horizontal plane and
%    the positive Z direction is upwards.  AZIMUTH angles are specified
%    relative to the positive X-axis [1 0 0], with positive azimuth angles
%    corresponding to positive Y-axis deviations.

%  jgd 4/11/97

% Assign default value, if necessary
	if nargin < 4,
		reference = [0 0 0];
	elseif nargin < 3,
		error('Not enough input arguments');
	end

%
% Orient input vectors and initialize output
%
	azimuth = azimuth(:);
	elev = elev(:);
	dist = dist(:);
	reference = ones(size(dist))*reference(:).';
	coord = zeros(length(azimuth),3);

% Get Z coord first
	coord(:,3) = dist.*sin(elev*pi/180);

% Get X and Y coords.
	dist_xy = dist.*cos(elev*pi/180);
	coord(:,1) = dist_xy.*cos(azimuth*pi/180);
	coord(:,2) = dist_xy.*sin(azimuth*pi/180);

% Make relative to reference;
	coord = coord + reference;