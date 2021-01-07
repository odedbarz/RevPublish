function h = shrinkage_fun(mean_strf, se_strf, gamma, shrinkage_type)
%
% function h = shrinkage_fun(mean_strf, se_strf, gamma, shrinkage_type, verbose)
%
%


if 4 > nargin
    shrinkage_type = 1;
end

% Choose the shrinkage function to use
switch shrinkage_type
    case 1
        % Stephen David et al., 2004, Natural Stimulus Statistics Alter the 
        % Receptive Field structure of V1 neurons
        snr = sqrt( max(0, 1 - gamma*(se_strf./mean_strf).^2) );
        h = mean_strf .* snr;

    case 2
        % *** Taken from STRFLAB, df_fast_filter_filter.m ***
        epsilon = 10^-8; % To prevent division by 0.
        snr = (1 + tanh(2.*(abs(mean_strf)-abs(gamma*se_strf))./(epsilon + abs(se_strf))))/2;
        h = snr .* mean_strf;

    otherwise
        error('ERROR in [STRF_C/shrinkage_fun.m]: unrecognized shrinkage_type!!!');
end






