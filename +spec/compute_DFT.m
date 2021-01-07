function X = compute_DFT(x, f, fs, dim)
%
% This function was taken from MAT:AB's signal processing toolbox, 
%   C:\Program Files\MATLAB\R2017a\toolbox\signal\signal\private\
% 
% I use it specifically for computinf DFT for nonuniform spaced
% frequencies, e.g., log scale frequencies.
%
%

if 4 > nargin
    dim = 1;
end

% see if we can get a uniform spacing of the freq vector
[~, ~, ~, maxerr] = spec.get_uniform_approx(f);

% see if the ratio of the maximum absolute deviation relative to the
% largest absolute in the frequency vector is less than a few eps
isuniform = maxerr < 3*eps(class(f));

if isuniform && size(x,1) == length(f)
    X = fft(x, [], dim);
    
    % Removes the conjugate spectrum part of the DFT Keep only the desired frequencies
    lowfreq = min(f);
    highfreq= max(f);

    f_full = linspace(0, fs, length(f))';
        
    f0_idx = lowfreq  <= f_full;
    f1_idx = highfreq >= f_full;
    X    = X(f0_idx & f1_idx, :);
        
else
    % for nonuniform (exponential) frequencies
    X = spec.compute_DFT_via_Goertzel(x, f, fs);
    
end

