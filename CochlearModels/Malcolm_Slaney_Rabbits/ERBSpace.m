function cfArray = ERBSpace(lowFreq, highFreq, N, rabbit_flag)
% function cfArray = ERBSpace(lowFreq, highFreq, N)
% This function computes an array of N frequencies uniformly spaced between
% highFreq and lowFreq on an ERB scale.  N is set to 100 if not specified.
%
% See also linspace, logspace, MakeERBCoeffs, MakeERBFilters.
%
% For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
% "Suggested formulae for calculating auditory-filter bandwidths and
% excitation patterns," J. Acoust. Soc. Am. 74, 750-753.
%
%
%
% 08/30-2018, Barzelay: added rabbit's cochlea
%

if nargin < 1
    % 08/30-2018, Barzelay: added rabbit's cochlea
	%lowFreq = 100;  
    lowFreq = 360;  
end

if nargin < 2
	highFreq = 44100/4;
end

if nargin < 3
	N = 100;
end


% 08/30-2018, Barzelay: added rabbit's cochlea
if nargin < 4
    rabbit_flag = false;
end

if ~rabbit_flag
    % Change the following three parameters if you wish to use a different
    % ERB scale.  Must change in MakeERBCoeffs too.
    EarQ = 9.26449;				%  Glasberg and Moore Parameters
    minBW = 24.7;
    %order = 1;

else    % *** Use Rabbit's ERB ***
    Fmin_rabbit = 360;    	% [Hz] minimum frequency of hearing range of rabbits    
    Fmax_rabbit = 42e3;   	% [Hz] maximum frequency of hearing range of rabbits    
    assert(lowFreq >= Fmin_rabbit, '--> Error at [ERBSpace.m]: lowFreq MUST BE greater or equal than %g Hz', Fmin_rabbit);
    assert(highFreq <= Fmax_rabbit, '--> Error at [ERBSpace.m]: highFreq MUST BE lower or equal than %g Hz', Fmax_rabbit);
    minBW = rabbitERB(Fmin_rabbit);
    f_lin = linspace(highFreq, lowFreq, N);     % desired frequencies (linear scale)
    [~, EarQ] = rabbitERB(f_lin);
    EarQ = EarQ(:);
end


% All of the followFreqing expressions are derived in Apple TR #35, "An
% Efficient Implementation of the Patterson-Holdsworth Cochlear
% Filter Bank."  See pages 33-34.
cfArray = -(EarQ*minBW) + exp((1:N)'.*(-log(highFreq + EarQ*minBW) + ...
        log(lowFreq + EarQ*minBW))/N) .* (highFreq + EarQ*minBW);








