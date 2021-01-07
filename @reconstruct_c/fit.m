function G = fit(obj, X_train, r_train, varargin)
%
%   function fit(obj, X_train, r_train)
%


%% Set the input data
p = inputParser;

% Required
% addRequired(p, 'obj', @isobject);            
addRequired(p, 'X_train', @isnumeric);            
addRequired(p, 'r_train', @isnumeric);            

% Optional
addOptional(p, 'lags_ms', obj.lags_ms, @isnumeric);         % (ms)
addOptional(p, 'f', obj.f, @isnumeric);                     % (Hz) frequency axis     
addOptional(p, 'algo_type', obj.algo_type, @isstr);           
addOptional(p, 'jk_flag', obj.jk_flag, @(x) isnumeric(x) || islogical(x));           
addOptional(p, 'n_splits', 11, @isnumeric);           
addOptional(p, 'gamma_4_JK', obj.ridge.gamma_4_JK, @isnumeric);           
addOptional(p, 'iscausal', obj.iscausal, @isnumeric);           
addOptional(p, 'inv_type', obj.ridge.inv_type, @isnumeric);           
addOptional(p, 'xcorr_type', obj.ridge.xcorr_type, @isnumeric);           
addOptional(p, 'verbose', 0, @isnumeric);           
addOptional(p, 'fignum', 0, @isnumeric);    

parse(p, X_train, r_train, varargin{:});
pars = p.Results;

% Update the jackknife flag, if needed
obj.f               = pars.f;
obj.lags_ms         = pars.lags_ms;     
obj.iscausal        = pars.iscausal;
obj.algo_type       = pars.algo_type;
obj.jk_flag         = pars.jk_flag;
obj.n_splits        = pars.n_splits;     
obj.ridge.inv_type  = pars.inv_type;
obj.ridge.xcorr_type= pars.xcorr_type;



%% 
if 0 == obj.jk_flag
    obj.X_train = X_train;
    obj.r_train = r_train;

    
    % Calculates 1 batch of filters
    obj.calc_filter(X_train, r_train,...
        'lags_ms', obj.lags_ms,...
        'f', obj.f,...        
        'algo_type', obj.algo_type,...        
        'iscausal', obj.iscausal,...
        'verbose', pars.verbose,...
        'fignum', pars.fignum );
    
    
else
    assert(~strcmpi('asd', obj.algo_type), '--> ASD currently does not work with JK!');
    
    % Apply jackknife to get the filters
    obj.jackknife(X_train, r_train,...
        'lags_ms', obj.lags_ms,...
        'n_splits', obj.n_splits,...
        'algo_type', obj.algo_type,...        
        'gamma_4_JK', obj.ridge.gamma_4_JK,...
        'iscausal', obj.iscausal,...
        'verbose', pars.verbose,...
        'fignum', pars.fignum );
    
end


%% Output
if 1 <= nargout
    G = obj.G;
    
end





