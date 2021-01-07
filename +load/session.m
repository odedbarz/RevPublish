function [S, n_sync, spikecount, fn_session] = session(input_par, verbose)
%
%   function [S, n_sync, spikecount, fn_session] = load.session(input_par, verbose)
%
% Input:
%    input_par: (struct) a measurement structure that contains the following
%    fiellds.
%         input_par.session   = 17;   
%         input_par.unit      = 1;   
%         input_par.measNum   = 5;    
%         input_par.syncchan  = 0;
%         input_par.spikechan = 1; 
%
% Output:
%   S         : (struct) Impale's structure
%   n_sync    : 
%   spikecount:
%
% Description:
%   Loads Impale's recorded structure.
%

import load.*

if 2 > nargin
    verbose = false;
end

% In case the input is is a structure
if isstruct(input_par)
    input_par = struct2table(input_par);
    assert(1 == size(input_par,1), '--> Currently, the input must contain ONE session to load!');
end

% In case the input is a table
if istable(input_par) 
    isfield_in_table = @(str) any(ismember(input_par.Properties.VariableNames, str));

    if isfield_in_table('syncchan')
        syncchan = input_par.syncchan;
    else
        syncchan = 0;
    end

    if isfield_in_table('spikechan')
        spikechan = input_par.spikechan;
    else
        spikechan = 1;
    end

    if iscell(input_par.ID)
        input_par.ID = input_par.ID{1};
    end

    if iscell(input_par.measType)
        input_par.measType = input_par.measType{1};
    end


    session_str = sprintf('OBZ_%s_%03d', input_par.ID, input_par.session);
    path2data = [path_to_data, session_str, filesep];
    fn_session = sprintf('OBZ_%s_%03d-%d-%d-%s_MERGED',...
        input_par.ID, input_par.session, input_par.unit, input_par.measNum, input_par.measType);
    fn_session = [path2data, fn_session, '.mat'];

    if verbose
        fprintf('--> Starting load.session...\n');
        fprintf('--> spikechan: %d\n', spikechan);
        fprintf('--> syncchan : %d\n', syncchan);
        fprintf('--> ID       : %d\n', input_par.ID);
        fprintf('--> measType : %d\n', input_par.measType);
    end


elseif ischar(input_par)
    assert( isempty(dir(input_par)), '--> Can''t find this file!' );
    
    % This is a valide (full) filename
    fn_session = input_par;
else
    error('--> Unrecognized input!!!');

end




%%

% The loaded session is already "re-ordered"
S = load.Impale_struct( fn_session );

% Extract the selected spikechan
spikes = cellfun(@(CH) CH==spikechan, S.ch, 'UniformOutput', false);
syncs  = cellfun(@(CH) CH==syncchan, S.ch, 'UniformOutput', false);

% Keep the desired spike-channel and sync-channel
t  = cellfun(@(T,IDX1,IDX2) T(IDX1 | IDX2), S.t, spikes, syncs, 'UniformOutput', false);
ch = cellfun(@(CH,IDX1,IDX2) CH(IDX1 | IDX2), S.ch, spikes, syncs, 'UniformOutput', false);

% Spike times
[tps, ~] = cellfun(@(T,CH) spiketools.pstimes(T,CH), t, ch, 'UniformOutput', false);


%% Keep the sync that corresponds to SPIKECHAN
% Only syncchan indices
syncs_slc_chan = cellfun(@(CH) CH == syncchan, ch, 'UniformOutput', false);   

% Only spikechan indices
spikes_slc_chan = cellfun(@(CH) CH == spikechan, ch, 'UniformOutput', false);

% Remove extra SYNCCHAN
spike_syncs = cellfun(@(SYNC,SPKS) [0; SYNC(1:end-1)] & SPKS, syncs_slc_chan, spikes_slc_chan, 'UniformOutput', false);
spike_syncs = cellfun(@(SPKS) [SPKS(2:end); false], spike_syncs, 'UniformOutput', false);

% The next loop will left the unused syncchan as -1
for kk = 1:length(ch(:))
    ch{kk}(syncs_slc_chan{kk}) = ch{kk}(syncs_slc_chan{kk}) -1;
    ch{kk}(spike_syncs{kk}) = ch{kk}(spike_syncs{kk}) +1;
end




%%
% Overwrite the selected SPIKECHAN
S.t   = t;
S.ch  = ch;
S.tps = tps;

n_sync = cellfun(@(CH) sum(0 == CH), S.ch, 'UniformOutput', true);
spikecount = cellfun(@(CH, T) sum(0 ~= CH), S.ch, 'UniformOutput', true);








