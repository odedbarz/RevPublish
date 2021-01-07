classdef strf_c < handle
    
    properties
        % Initialize the output variables
        binwidth    = []       % (ms)
        f           = []       % (Hz) frequency axis used for plotting ONLY
        X_train     = []
        r_train     = []        
        r_est       = []  
        
        theta       = []       % shrinkage coefficient
        
        % Ridge regression
        gamma       = 0.0005   % tolerance in choosing the eigenvalues
        
        % jackknife
        jk_flag     = 0
        theta_4_JK  = 1:0.25:2       % list of shrinkage options
        
        n_splits    = 11        
        bias_choice = 1        % remove the bias
        
        % xcorr_stim_psth files these fields
        bias        = struct('Sft', [], 'psth', [], 'bias', [])   
        
        Ravg        = []
        Rw          = []        
        Cavg        = []
        Cw          = []
        
        % Algorithm type to calculate the STRF 
        algo_type   = 'regression';     % {'ASD', 'regression'}
        asdstats    = []
        ridge       = []        % 1: use convmtx, 2: use for-loop (Nima's lab) procedures        
        xcorr_type  = 3         % (1x1) correlation methos
        strf        = []
        strf_noncausal          % save the whole response too. If iscasual == 1, 
                                % then strf equals strf_noncausal
        
        % Add a window to smooth the response(especially for the PSTH)
        smooth_response_flag = 0
        
        % Computation parameters (default values)
        gpu_flag    = 0 < gpuDeviceCount   % use the GPU, if available
        
        %verbose     = 0
        %fignum      = []
    end
    
    properties (SetObservable)
        lags_ms   = 70         % (ms)
        iscausal = false       % if true, the reconstructed filters are bewteen 0 ms and lag_ms         
    end
    
    properties (SetAccess = private)
        lags   = []
        n_lags = []
    end
    
    
    methods
            % *******************
            % *** CONSTRUCTOR ***
            % *******************            
            function obj = strf_c( binwidth, varargin )
                
                p = inputParser;
                
                % *** Required parameters ***
                % This is the x-axis time resolution
                addRequired(p, 'binwidth', @isnumeric);             % (ms)        

                
                % *** Optional parameters ***
                addOptional(p, 'f', [], @isnumeric);                % (Hz)
                
                % # of lags for xcorr. N_WIN must be even. if given as an odd number, it will 
                %   be reduced by one to an even number.                
                addOptional(p, 'lags_ms', obj.lags_ms, @isnumeric); 
                addOptional(p, 'iscausal', obj.iscausal, @(x) islogical(x) || isscalar(x)); 
                
                % Algorithm type to calculate the STRF 
                addOptional(p, 'algo_type', obj.algo_type, @isstr);      % % {'ASD', 'regression'}             
                addOptional(p, 'jk_flag', obj.jk_flag, @isnumeric);  
                addOptional(p, 'n_splits', obj.n_splits, @isnumeric); 
                
                % Tolerance for the ridge regression. If tol is an array, the function will 
                %   return an array (s) of STRFs that relates to each of the M values.                
                addOptional(p, 'gamma', 0.0005, @isnumeric);   
                addOptional(p, 'gamma_4_JK', [], @isnumeric);   
                                
                addOptional(p, 'bias_choice', obj.bias_choice, @isnumeric);  
                addOptional(p, 'xcorr_type', obj.xcorr_type, @isnumeric);  
                addOptional(p, 'gpu_flag', obj.gpu_flag, @isnumeric);  

                % For the plottings ONLY
                addOptional(p, 'faxis', [], @isnumeric);        % (n_bands x 1) frequency axis (Hz)   
                addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
                %addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

                parse(p, binwidth, varargin{:});
                pars = p.Results;
                                
                % Update all input values into the object
                obj.binwidth    = binwidth;             % (ms)
                obj.f           = pars.f;               % (Hz)
                obj.iscausal    = pars.iscausal;
                obj.xcorr_type  = pars.xcorr_type;
                obj.jk_flag     = pars.jk_flag;    
                obj.n_splits    = pars.n_splits;    
                obj.algo_type   = pars.algo_type;
                obj.gpu_flag    = pars.gpu_flag;

                % Create & update the ridge structure
                obj.init_ridge_regression_pars(pars);
                
                addlistener(obj, 'lags_ms', 'PostSet', @obj.update_lags);
                addlistener(obj, 'iscausal', 'PostSet', @obj.update_lags);
                obj.lags_ms = pars.lags_ms;                     
                
                if pars.verbose
                    fprintf('--> binwidth   : %d ms\n', obj.binwidth);
                    fprintf('--> lags_ms    : %g ms\n', obj.lags_ms);
                    fprintf('--> inv_type   : [%s]\n',  num2str(obj.inv_type, ' %g'));
                    fprintf('--> bias_choice: %d\n', obj.bias_choice);
                end
            end
            
            % Initialize ridge regression structure
            function init_ridge_regression_pars(obj, pars)
                
                if strcmpi('ASD', obj.algo_type)
                    % If algorithm is ASD, keep the ridge field empty                    
                    obj.ridge = struct('gamma', [], 'gamma_4_JK', [], 'xcorr_type', []);
                    return;
                end
                
                % default value; JK finds a better value using crossover
                if isfield(pars, 'gamma') && ~isempty(pars.gamma)
                    st.gamma = pars.gamma;	% user input 
                else
                    st.gamma = 0.0005;   	% (default)
                end
                
                % TOLERANCE list for the ridge regression (JK)
                if isfield(pars, 'gamma_4_JK') && ~isempty(pars.gamma_4_JK)
                    st.gamma_4_JK = pars.gamma_4_JK;	% user input 
                else
                    st.gamma_4_JK = logspace(log10(1e-7), log10(0.1), 15);   	% (default)
                end
                
                % Apdate the object
                obj.ridge = st;
            end
                        
    end
    
    
    % DEPENDENT properties
    properties (Dependent)
        n_bands     % # of frequency bands
        nt          % # samples along the time domain
        fs          % (Hz) sampling rate along the time-axis
        n_units     % # of units used 
        bf          % (1x1) best frequency of the STRF (maximum frequency)
    end
    
    
    % DEPENDENT methods       
    methods
        % # of frequency bands
        function n = get.n_bands(obj)
            n = size(obj.X_train, 1);
        end
        
        % # samples along the time domain
        function n = get.nt(obj)
            n = size(obj.X_train, 2);
        end
        
        % (Hz) sampling rate along the time-axis
        function Fs = get.fs(obj)
            Fs = 1/(1e-3*obj.binwidth);
        end
        
        function n = get.n_units(obj)
            n = size(obj.r_train, 2);
        end
       
        function bf = get.bf(obj)
            if ~isempty(obj.strf)
                bf = obj.calc_bf;
            else
                bf = nan;
            end
        end
        
    end
    
    
    % Event Methods
    methods 
        function update_lags(obj, src, evnt)
            %
            %   update_lags(obj, src, evnt)
            %
            % Description: 
            % Perform this every time that the LAGS_MS or ISCASUAL is updated.
            one_half_smp   = ceil(obj.lags_ms / obj.binwidth); 
            
            % Make sure that N_WIN is even
            one_half_smp = one_half_smp - mod(one_half_smp,2);
            
            % Lag vectos for the reconstructed matrix
            if obj.iscausal        
                lags_v = 0:one_half_smp;   
            else
                lags_v = (-one_half_smp):one_half_smp;                
            end
            assert(lags_v(1)<=0, '--> Check out that the assigned LAGS are valid!');
           
            
            % Update the object
            obj.lags = lags_v;
            obj.n_lags = numel(lags_v);            
        end
    end
    

    
    
end

