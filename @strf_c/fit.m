function [strf, JN] = fit(obj, X_train, r_train, varargin)
%
%   function fit(obj, X_train, r_train)
%


%% Set the input data
p = inputParser;

% Required
% addRequired(p, 'obj', @isobject);            
addRequired(p, 'X_train', @isnumeric);            
addRequired(p, 'r_train', @isnumeric);            

% Optional -- single STRF
addOptional(p, 'iscausal', obj.iscausal, @isnumeric);           
addOptional(p, 'gamma', obj.ridge.gamma, @isnumeric);      % (1x1) tolerance for the ridge regression
addOptional(p, 'xcorr_type', obj.xcorr_type, @isnumeric);    % (1x1)  cross-correlation type to use

% Algorithm type to calculate the STRF 
addOptional(p, 'algo_type', obj.algo_type, @istr);      % % {'ASD', 'regression'}             

% Optional -- jackknife
addOptional(p, 'jk_flag', obj.jk_flag, @(x) isnumeric(x) || islogical(x));           
addOptional(p, 'n_splits', 11, @isnumeric);           
addOptional(p, 'gamma_4_JK', obj.ridge.gamma_4_JK, @isnumeric);           
addOptional(p, 'theta_4_JK', obj.theta_4_JK, @isnumeric);           

% Optional -- plotting & debug
addOptional(p, 'verbose', 0, @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

parse(p, X_train, r_train, varargin{:});
pars            = p.Results;


% Update the jackknife flag, if needed
obj.jk_flag    = pars.jk_flag;
obj.algo_type  = pars.algo_type;
obj.xcorr_type = pars.xcorr_type;
obj.ridge.gamma_4_JK = pars.gamma_4_JK;
obj.theta_4_JK = pars.theta_4_JK;
obj.iscausal   = pars.iscausal;
obj.n_splits   = pars.n_splits;     



%% 
if 0 == pars.jk_flag
    % Apply ridge regression
    strf = obj.calc_strf(X_train, r_train,...
        'gamma', obj.gamma, ...
        'xcorr_type', obj.xcorr_type, ...
        'fignum', pars.fignum, 'verbose', pars.verbose);
    
    JN = [];
    
else
    % Apply jackknife to get the filters
    [strf, JN] = obj.jackknife(X_train, r_train,...
        'xcorr_type', obj.xcorr_type, ...     
        'iscausal', obj.iscausal,...        
        'n_splits', obj.n_splits,...
        'gamma_4_JK', obj.ridge.gamma_4_JK,...
        'theta_4_JK', obj.theta_4_JK,...
        'fignum', pars.fignum, 'verbose', pars.verbose);
    
    
end






