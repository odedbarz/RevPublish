function [psth_mtx, bins, data_st] = calc_PSTHs(S_list, binwidth, varargin)
%
% function [psth_mtx, bins, data_st] = calc_PSTHs(S_list, binwidth, varargin)
%
% Inputs:
%   S_list
%   binwidth
%   [pst_type]
%   [syncchan]
%   [n_stimuli]
%   [verbose]
%% Set the inputs
p = inputParser;


% (table) A table that holds all relevant information about the neurons in
%         the database.
addRequired(p, 'S_list', @(x) iscell(x) || isstruct(x));     
addRequired(p, 'binwidth', @isnumeric);   % (ms)         

addOptional(p, 'pst_type', {'pstw', 50}, @iscell);      % {'psth'}, {'pstw', WIN_SIZE_MS}
addOptional(p, 'syncchan', 0, @isnumeric);    
% addOptional(p, 'n_stimuli', 5, @isnumeric);         % # of DRR stimuli (# of columns(second dim.) in the PSTH matrix
addOptional(p, 'verbose', [], @isnumeric);          % (1x1) write things to the command line?
% addOptional(p, 'fignum', [], @isnumeric);         % (1x1) if not empty, plot stuff (for debug)

parse(p, S_list, binwidth, varargin{:});


pst_type  = p.Results.pst_type;       
syncchan  = p.Results.syncchan;       
% n_stimuli = p.Results.n_stimuli;       
% verbose     = p.Results.verbose;       
% fignum   = p.Results.fignum;       

drr = get_DRR_list_and_indices;
n_stimuli = length(drr.labels);


if isstruct(S_list)
    S_list = {S_list};
end

data_st.binwidth = binwidth;
data_st.pst_type = pst_type;
data_st.syncchan = syncchan;


%% Init. analysis table
n_neurons   = length(S_list);
duration_ms = S_list{1}.info.duration_ms;   % all PSTHs must have the same duration!
n_smp       = duration_ms/binwidth;
psth_mtx    = nan(n_smp, n_stimuli, n_neurons);
psth_list   = cell(1, n_neurons);


for k = 1:n_neurons    
    S = S_list{k};

    assert(duration_ms == S.info.duration_ms);
    spikechan = S.info.spikechan;
    
    %% ** PSTH **     
    switch upper(pst_type{1})
        case 'PSTH'
            psth_k = PSTH(S,...
                duration_ms,...
                binwidth,...   % !! set the PSTH bins to equal that of the spectrogram
                spikechan,...
                syncchan );
            
        case 'PSTW'    
            pstw_win_size_ms = pst_type{2};     
            psth_k = PSTW(S,...
                duration_ms,...
                binwidth,...   % !! set the PSTH bins to equal that of the spectrogram
                pstw_win_size_ms, ...
                spikechan,...
                syncchan );
    
        otherwise
            error('--> Error in [load.response.m]: no such PST function!');
    end
    
%     %assert( nnz(S.info.n_sync) == n_stimuli );
%     if nnz(S.info.n_sync) ~= n_stimuli
%         fprintf('--> ** skipping the %d measurement!! **\n', k);
%         continue;
%     end
    
    % 
    psth_mtx(:, psth_k.available_meas_idx, k) = psth_k.H(:,psth_k.available_meas_idx);
    psth_list{k} = psth_k;
    
end

bins = psth_list{1}.bins;       % (ms)
data_st.psth_list = psth_list;















