function m = gen_array(M,array_type,args)
%
%	m = gen_array(M,array_type,args)
%
%	Function to return the (x,y,z) coordinates, based upon an array center of
%	(0,0,0), for an M-element array described by the parameters in args.
%	All arrays considered here are two-dimensional and corrdinate positions
%	are given for the horizontal (XY) plane only.
%
%	Array types include:
%
%		'circ1'		M-elements equally-spaced along the rim of a
%				circle of radius r meters.
%				args = [r];
%
%		'circ2'		M-elements divided into two sets of M/2 elements.
%				Each set of M/2 elements is equally spaced along
%				the rim of a circle -- with radii of r1 and r2 
%				meters, respectively.
%				args = [r1; r2];
%
%		'circ3		M-elements divided into two sets of M/2 elements.
%				Both sets of elements are equally-spaced (within 
%				each set) along the same circle of radius r 
%				meters.  The second set of array elements is 
%				rotated about the center by an angle of th degrees.
%				args = [r; th];
%
%		'cross1'	M-elements divided into two sets of M/2 elements.
%				Each set of M/2 element is arranged into a ULA
%				architecture with a spacing of d meters.  The
%				two ULAs are arranged into a perfect cross --
%				90-degree inter array angles, crossing at ULA 
%				centers.
%				args = [d];
%
%		'cross2'	M-elements divided into two sets of M/2 elements.
%				Each set of M/2 element is arranged into a 
%				linear array with spacings that alternate between
%				d1 and d2 meters.  The two linear arrays are 
%				arranged into a perfect cross -- 90-degree inter 
%				array angles, crossing at array centers.
%				args = [d1; d2];
%
%		'ULA' M-elements arranged as a uniformly-spaced linear
%				array with inter-element spacing of d meters.
%				Array is centered at the origin and makes an angle of
%				th degrees with the positive x axis.
%				args = [d; th]
%


	if strcmp(array_type,'circ1'),

		% Circular array Type I
		r = args(1);
		th_step = pi/(M/2);
		m = r*[cos(th_step*[0:M-1]'), sin(th_step*[0:M-1]'), ...
			zeros(M,1)];
		
	elseif strcmp(array_type,'circ2'),

		% Circular array Type II
		r1 = args(1); 
		r2 = args(2);
		th_step = pi/(M/4);
		m = r1*[cos(th_step*[0:M/2-1]'), sin(th_step*[0:M/2-1]'), ...
			zeros(M/2,1)];
		m = [m;r2*[cos(th_step*[0:M/2-1]'), sin(th_step*[0:M/2-1]'), ...
			zeros(M/2,1)]];

	elseif strcmp(array_type,'circ3'),

		% Circular array Type III
		r = args(1); 
		th = args(2)*pi/180;
		th_step = pi/(M/4);
		m = r*[cos(th_step*[0:M/2-1]'), sin(th_step*[0:M/2-1]'), ...
			zeros(M/2,1)];
		m = [m;r*[cos(th_step*[0:M/2-1]'+th), sin(th_step*[0:M/2-1]'+th), ...
			zeros(M/2,1)]];

	elseif strcmp(array_type,'cross1'),

		% Cross array Type I
		d = args(1);
		span = (M/2-1)*d;
		m = [[zeros(M/2,1),[-span/2:span/(M/2-1):span/2]',zeros(M/2,1)]; ...
			[[-span/2:span/(M/2-1):span/2]',zeros(M/2,1),zeros(M/2,1)]];

	elseif strcmp(array_type,'cross2'),

		% Cross array Type II
		d1 = args(1);
		d2 = args(2);
		d = zeros(M/2-1,1); 
		d(1:2:length(d)) = d1;
		d(2:2:length(d)) = d2;
		d = cumsum([0;d]);
		span = d1*floor(M/4) + d2*ceil(M/4-1);
		m = [[zeros(M/2,1),-span/2+d,zeros(M/2,1)]; ...
			[-span/2+d,zeros(M/2,1),zeros(M/2,1)]];
   
	elseif strcmp(array_type,'ULA'),
   
   	% Uniform Linear Array
      d = args(1);
		th = args(2);
		span = (M-1)*d;
		m = [[-span/2:span/(M-1):span/2]',zeros(M,1),zeros(M,1)];
  		Rot = [cos(th*pi/180) sin(th*pi/180) 0; ...
				-sin(th*pi/180) cos(th*pi/180) 0; ...
				0 0 1];
		m = m*Rot;
	end
