function [rates_mtx, sem_mtx] = calc_Rates(S_list, t_start, t_end, varargin)
%
% function [rates_mtx, sem] = calc_Rates(S_list, binwidth, varargin)
%
% Inputs:
%   S_list
%   t_start, t_end
%   [isi_order]
%   [n_stimuli]
%   [syncchan]
%   [verbose]
%% Set the inputs
p = inputParser;


% (table) A table that holds all relevant information about the neurons in
%         the database.
addRequired(p, 'S_list', @(x) iscell(x) || isstruct(x));     
addRequired(p, 't_start', @isnumeric);              % (ms)   
addRequired(p, 't_end', @isnumeric);                % (ms)   

addOptional(p, 'n_stimuli', 5, @isnumeric);         % # of DRR stimuli (# of columns(second dim.) in the PSTH matrix
addOptional(p, 'syncchan', 0, @isnumeric);    
addOptional(p, 'verbose', [], @isnumeric);          % (1x1) write things to the command line?
% addOptional(p, 'fignum', [], @isnumeric);         % (1x1) if not empty, plot stuff (for debug)

parse(p, S_list, t_start, t_end, varargin{:});


syncchan    = p.Results.syncchan;       
n_stimuli   = p.Results.n_stimuli;       
% verbose     = p.Results.verbose;       
% fignum   = p.Results.fignum;       

if isstruct(S_list)
    S_list = {S_list};
end


%% Init. analysis table
n_neurons   = length(S_list);
duration_ms = S_list{1}.info.duration_ms;   % all PSTHs must have the same duration!
t_start     = 0;                % (ms)
t_end       = duration_ms;      % (ms)

rates_mtx   = zeros(n_neurons, n_stimuli);
sem_mtx     = zeros(n_neurons, n_stimuli);

for k = 1:n_neurons    
    S = S_list{k};

    % Make sure that all measurements have the same duration!
    assert(duration_ms == S.info.duration_ms);
    
    [rates_k, sem_k] = Rates(t_start, t_end,...
                        S_list{k}.t, S_list{k}.ch,...
                        'spikechan', S.info.spikechan,...
                        'syncchan', syncchan );

    rates_mtx(k,:) = rates_k;
    sem_mtx(k,:) = sem_k;
    
end



