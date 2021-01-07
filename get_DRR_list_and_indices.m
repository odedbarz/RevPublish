function drr = get_DRR_list_and_indices(verbose)
%
%   function drr = get_DRR_list_and_indices(verbose)
%
% * The ordering of the stimuli in the Impale files
%     Revb.            20%        80%   100%        20%        80%     100%           
%     Dist.           1.5m       1.5m   1.5m       3.0m       3.0m     3.0m
%
% * The NEW ordering is:
%     #         1        2       3        4        5       6
%     DRR      Dry  9.4 dB  4.8 dB  -2.5 dB  -8.2 dB      --
%     Revb.   100%     80%     80%      20%      20%    100%           
%     Dist.   1.5m    1.5m    3.0m     1.5m     3.0m    3.0m

%

if 0 == nargin
    verbose = false;
end

% Revb.            20%        80%   100%        20%        80%     100%           
% Dist.           1.5m       1.5m   1.5m       3.0m       3.0m     3.0m
drr.labels= {'-2.5 dB',  '9.4 dB', 'Dry', '-8.2 dB',  '4.8 dB',  'Dry'};
drr.dist  = [      1.5,       1.5,    1.5,       3.0,       3.0,    3.0];
drr.revb  = [       20,        80,    100,        20,        80,    100];
drr.dry    = 3;
drr.dry2   = 6;     % DRY condition, second option
drr.sortby = [drr.dry, 2, 5, 1, 4, 6];     % ordered from high (min reverb) to low (high reverb) DRR
%drr.labels = drr_list(drr.sortby);         % labels are ordered
if verbose
    fprintf('\n Labels:\n');
    disp(drr);
end

drr.n_drr = 5; %length(drr.labels);
drr.ordered = drr.sortby(1:5);