function h = shrinkage_fun(X_mean, X_sd, gamma, shrinkage_type)
%
% function h = shrinkage_fun(X_mean, X_sd, gamma, shrinkage_type, verbose)
%
%


if 4 > nargin
    shrinkage_type = 1;
end

% Choose the shrinkage function to use
switch shrinkage_type
    case 1
        %   David, S.V., Vinje, W.E. and Gallant, J.L., 2004. Natural stimulus 
        %   statistics alter the receptive field structure of v1 neurons. Journal of 
        %   Neuroscience, 24(31), pp.6991-7006.
        snr = sqrt( max(0, 1 - gamma*(X_sd./X_mean).^2) );
        h = X_mean .* snr;

    case 2
        % *** Taken from STRFLAB, df_fast_filter_filter.m ***
        epsilon = 10^-8; % To prevent division by 0.
        snr = (1 + tanh(2.*(abs(X_mean)-abs(gamma*X_sd))./(epsilon + abs(X_sd))))/2;
        h = snr .* X_mean;

    otherwise
        error('ERROR in [strfpkg.shrinkage_fun.m]: unrecognized shrinkage_type!!!');
end






