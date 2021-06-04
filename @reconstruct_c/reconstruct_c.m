classdef reconstruct_c < handle
    
    properties
        % Initialize the output variables
        binwidth    = [];           % (ms)
        f           = [];           % (Hz) frequency axis; used for plotting ONLY
        X_train     = [];
        r_train     = [];        
        X_est       = [];        
        jk_flag     = 0;            % jackknife flag
        n_splits    = 11;        
        eigv        = [];
        %bias_choice = 1;            % remove the bias
        %bias        = [];           % (struct) hold the spectrogram & responses biases
        %inv_type    = 2;    %[2, 0.005];
        G           = [];
        G_sd        = [];
        Crr         = [];
        Crr_inv     = [];
        Crs         = [];
        R_test      = [];
        ridge       = struct('gamma', [], 'gamma_4_JK', [],...
                        'xcorr_type', [], 'inv_type', []);
        bias        = struct('remove', 1, 'train', [], 'test', []);
        
        
        % Computation parameters (default values)
        %algo_type   = 1     % 1: use convmtx, 2: use for-loop (Nima's lab) procedures
        gpu_flag    = 0     %0 < gpuDeviceCount;   % use the GPU, if available
        
        % Algorithm type to calculate the STRF 
        algo_type   = 'regression';     % {'ASD', 'regression'}        
        asdstats    = [];   % if used, hold statistics about the ASD
        
        %verbose     = 0;
        %fignum      = [];
    end
    
    
    properties (SetObservable)
        lags_ms   = 70;           % (ms)
        iscausal = false;        % if true, the reconstructed filters are bewteen 0 ms and lag_ms;         
    end
    
    properties (SetAccess = private)
        lags   = [];
        n_lags = [];
    end
    
    
    methods
            % *******************
            % *** CONSTRUCTOR ***
            % *******************            
            function obj = reconstruct_c( binwidth, varargin ) 
                p = inputParser;
                addRequired(p, 'binwidth', @isnumeric);             % (ms)        
                
                addOptional(p, 'f', [], @isnumeric);                % (Hz)
                addOptional(p, 'lags_ms', obj.lags_ms, @isnumeric); 
                addOptional(p, 'iscausal', obj.iscausal, @(x) islogical(x) || isscalar(x)); 
                addOptional(p, 'inv_type', [], @isnumeric);       
                addOptional(p, 'algo_type', obj.algo_type, @isstr);       
                addOptional(p, 'gamma', 0.005, @isnumeric);       
                addOptional(p, 'gamma_4_JK', logspace(log10(1e-5), log10(1.0), 14), @isnumeric);       
                addOptional(p, 'remove_bias', obj.bias.remove, @isnumeric);  
                addOptional(p, 'xcorr_type', 1, @isnumeric);  
                addOptional(p, 'gpu_flag', obj.gpu_flag, @isnumeric);  
                
                addOptional(p, 'verbose', 0, @isnumeric);           
                
                parse(p, binwidth, varargin{:});
                pars = p.Results;
                
                obj.binwidth    = binwidth;             % (ms)
                obj.f           = pars.f;               % (Hz)
                obj.iscausal    = pars.iscausal;
                obj.bias.remove = pars.remove_bias;  
                obj.gpu_flag    = pars.gpu_flag;                  

                % Algorithm type to calculate the STRF 
                obj.algo_type  = pars.algo_type;     % {'ASD', 'regression'}

                
                % Create & update the ridge structure
                obj.init_ridge_regression_pars(pars);
                
                addlistener(obj, 'lags_ms', 'PostSet', @obj.update_lags);
                addlistener(obj, 'iscausal', 'PostSet', @obj.update_lags);
                obj.lags_ms         = pars.lags_ms;                     
                
                if pars.verbose
                    fprintf('--> binwidth   : %d ms\n', obj.binwidth);
                    fprintf('--> lags_ms    : %g ms\n', obj.lags_ms);
                    fprintf('--> inv_type   : [%s]\n',  num2str(obj.ridge.inv_type, ' %g'));
                    fprintf('--> remove_bias: %d\n', obj.bias.remove);
                    fprintf('--> verbose    : %d\n', obj.verbose);
                    fprintf('--> fignum     : %d\n', obj.fignum);
                end
            end
            
            % Initialize ridge regression structure
            function init_ridge_regression_pars(obj, pars)

                if strcmpi('ASD', obj.algo_type)
                    % If algorithm is ASD, keep the ridge field empty                    
                    obj.ridge = struct('gamma', [], 'gamma_4_JK', [],...
                        'xcorr_type', [], 'inv_type', []);
                    return;
                end
                
                % default value; JK finds a better value using crossover
                if isfield(pars, 'gamma') && ~isempty(pars.gamma)
                    st.gamma = pars.gamma;	% user input 
                else
                    st.gamma = 0.005;   	% (default)
                end
                
                % GAMMA values for the JK
                if isfield(pars, 'gamma_4_JK') && ~isempty(pars.gamma_4_JK)
                    st.gamma_4_JK = pars.gamma_4_JK;	% user input 
                else
                    st.gamma_4_JK = logspace(log10(1e-5), log10(1.0), 14);   	% (default)
                end
                            
                % corr matrix algorithm 
                %   1: use convmtx
                %   2: use for-loop (Nima's lab) procedures
                %   3: another version that calculates the corr matrix by
                %      parts (to avoid memory constraints)
                if isfield(pars, 'xcorr_type') && ~isempty(pars.xcorr_type)
                    st.xcorr_type = pars.xcorr_type;	% user input 
                else
                    st.xcorr_type = 1;   	% (default)
                end
                
                
                % Default inversion type for the autocorrelation matrix
                if isfield(pars, 'inv_type') && ~isempty(pars.inv_type)
                    st.inv_type = pars.inv_type;	% user input 
                else
                    st.inv_type = 2;   	% (default)
                end
                
                % Apdate the object
                obj.ridge = st;
            end
            
            % Extract the n'th reconstruction filter
            function Gi = Gi(obj, n)
                if isempty(obj.G)
                    error('--> [@reconstruct_c->Gi]: There are no filters in this object (G is empty)!!!');
                elseif 0 >= n || n > obj.n_neurons
                    error('--> [@reconstruct_c->Gi]: the requested n''th (n = %d) filter isn''t available (try 1:%d)!!!', n, obj.n_neurons);
                end
                Gi = obj.G( (n-1)*obj.n_lags + (1:obj.n_lags), : );
            end
            
    end
    
    
    % DEPENDENT properties
    properties (Dependent)
        n_bands     % # of frequency bands
        nt          % # samples along the time domain
        n_neurons   % # of neurons\sessions
        fs          % (Hz) sampling rate along the time-axis
    end
    
    
    % DEPENDENT methods       
    methods
        % # of neurons participated in the fitting of the reconstruction filters
        function n = get.n_neurons(obj)
            if ~isempty(obj.G)
                n = size(obj.G,1)/obj.n_lags;
            elseif ~isempty(obj.r_train)
                n = size(obj.r_train, 2);
            else
                error('--> Class [reconstruction]: can''t find # of neurons!');
            end
            assert(fix(n)==n, '--> RECONSTRUCT_C.m: N_NEURONS must be an integer!!!');
            
        end
        
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
            
            % Lag vectos for the reconstructed matrix
            if obj.iscausal        
                lags_v = 0:one_half_smp;   
            else
                lags_v = -one_half_smp:one_half_smp;                
            end
            assert(lags_v(1)<=0, '--> Check out that the assigned LAGS are valid!');
            
            % Update the object
            obj.lags = lags_v;
            obj.n_lags = numel(lags_v);            
        end
    end
    

    
    
end

