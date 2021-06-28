function [isi_mtx, bins, data_st] = calc_ISI(S_list, binwidth, varargin)
%
% function [psth_mtx, data_st] = calc_PSTHs(S_list, binwidth, varargin)
%
% Inputs:
%   S_list
%   binwidth
%   [isi_order]
%   [syncchan]
%   [n_stimuli]
%   [verbose]
%% Set the inputs
p = inputParser;


% (table) A table that holds all relevant information about the neurons in
%         the database.
addRequired(p, 'S_list', @(x) iscell(x) || isstruct(x));     
addRequired(p, 'binwidth', @isnumeric);             % (ms)   
addOptional(p, 'n_bins', 100, @isnumeric); 
addOptional(p, 'isi_order', 1, @isnumeric);     
addOptional(p, 'n_stimuli', 5, @isnumeric);         % # of DRR stimuli (# of columns(second dim.) in the PSTH matrix
addOptional(p, 'syncchan', 0, @isnumeric);    
addOptional(p, 'verbose', [], @isnumeric);          % (1x1) write things to the command line?
% addOptional(p, 'fignum', [], @isnumeric);         % (1x1) if not empty, plot stuff (for debug)

parse(p, S_list, binwidth, varargin{:});


isi_order   = p.Results.isi_order;       
syncchan    = p.Results.syncchan;       
n_bins       = p.Results.n_bins;       
n_stimuli       = p.Results.n_stimuli;       
% verbose     = p.Results.verbose;       
% fignum   = p.Results.fignum;       

if isstruct(S_list)
    S_list = {S_list};
end

data_st.binwidth = binwidth;
data_st.pst_type = isi_order;
data_st.pst_type = n_bins;
data_st.pst_type = n_stimuli;
data_st.syncchan = syncchan;


%% Init. analysis table
n_neurons   = length(S_list);
duration_ms = S_list{1}.info.duration_ms;   % all PSTHs must have the same duration!
isi_mtx     = zeros(n_bins, n_stimuli, n_neurons);

for k = 1:n_neurons    
    S = S_list{k};

    % Make sure that all measurements have the same duration!
    assert(duration_ms == S.info.duration_ms);
    
    [isi_k, bins] = ISI(S_list{k}.t, S_list{k}.ch, 0, duration_ms,...
        'binwidth', binwidth,...
        'nbins', n_bins ,...
        'isi_order', isi_order,...
        'spikechan', S.info.spikechan,...
        'syncchan', syncchan ...
        );

    if size(isi_k,2) ~= n_stimuli
        continue;
    end
    
    isi_mtx(:,:,k) = isi_k;
    
end



