%	script CAL_DFR calculates
%	the diffuse-filed response for a sphere.

%	10-3-00...pmz


%	third-oct frequencies
	f = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 ...
           2500 3150 4000 5000 6300 8000];

%	f = (0:100:10000);

	c = 345;
	k = 2*pi*f/c;
	theta = (0:179)/180*pi;
	wt = sin(theta);	% sine weighting for spherical volume 

%	for 8.5 cm sphere
	rad = 0.085;
	yc_sph = sp(k,theta,rad,rad);
	pow_sph = yc_sph.*conj(yc_sph);
	pow_sph_wtd = diag(wt) * pow_sph;
	
%	for free field (very small radius sphere)
	rad = 0.000001;
	yc_ff = sp(k,theta,rad,rad);
	pow_ff = yc_ff.*conj(yc_ff);
	pow_ff_wtd = diag(wt) * pow_ff;

	dfr = 10*log10(sum(pow_sph_wtd)./sum(pow_ff_wtd));

	semilogx(f,dfr);