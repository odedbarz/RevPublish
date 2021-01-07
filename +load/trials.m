function [trials, S, tps, ch, fn2load] = trials(meas_st, verbose)
%
%   function [trial_st, S, fn2load] = load.trials(meas)
%
% Input:
%    meas_st: (struct) a measurement structure that contains the following
%    fiellds.
%         meas_st.session   = 17;   % 17; 17; 19
%         meas_st.unit      = 1;   %  1;  2;  2
%         meas_st.measNum   = 5;    %  5;  6; 10
%         meas_st.syncchan  = 0;
%         meas_st.spikechan = 1; 
% 
%         % Session indices to load:
%         meas_st.I = 1;
%         meas_st.J = 1;
%   
%   verbose: (logical) plots on screen.
%
% Output:
%   trial_st: (struct)
% 
%      spk_array: {1×N cell} N trials of spikes
%            tps: [K×1 double] K are number of spikes
%          nstim: N # of trials
%     trials_idx: [K×1 logical]
%              t: [K×1 double]
%             ch: [K×1 double]
%
% Description:
%   Loads spikes of a specified measurement. The spikes are ordered in
%   trials.
%

import load.*

if 2 > nargin
    verbose = false;
end

syncchan = meas_st.syncchan;
spikechan = meas_st.spikechan;

session_str = sprintf('OBZ_C74_%03d', meas_st.session);
path2data = ['../../.data/', session_str, '/'];
fn2load = sprintf('OBZ_C74_%03d-%d-%d-%s',...
    meas_st.session, meas_st.unit, meas_st.measNum, meas_st.measType);

% The loaded session is already "re-ordered"
S = load.Impale_struct( [path2data, fn2load, '.mat'] );
% assert( isfield(S.info, 'reordered') && S.info.reordered,...
%     '--> ERROR: The loaded session MUST be "re-ordered"!');


% Extract the selected spikechan
spikes = cellfun(@(CH) CH==spikechan, S.ch, 'UniformOutput', false);
chans  = cellfun(@(CH) CH==syncchan, S.ch, 'UniformOutput', false);

% Keep the desired spike-channel and sync-channel
t  = cellfun(@(T,IDX1,IDX2) T(IDX1 | IDX2), S.t, spikes, chans, 'UniformOutput', false);
ch = cellfun(@(CH,IDX1,IDX2) CH(IDX1 | IDX2), S.ch, spikes, chans, 'UniformOutput', false);

% Spike times
[tps, ~] = cellfun(@(T,CH) spiketools.pstimes(T,CH), t, ch, 'UniformOutput', false);


num_rep = cellfun(@(CH) nnz(syncchan == CH), ch, 'UniformOutput', true);
trials_idx  = cellfun(@(CH) cumsum( CH == syncchan ), ch, 'UniformOutput', false);
trials = cell(size(num_rep));

for ii = 1:size(num_rep,1)
    for jj = 1:size(num_rep,2)
        if isempty(tps{ii,jj}), continue; end
        trials{ii,jj} = arrayfun(@(IDX) tps{ii,jj}(IDX==trials_idx{ii,jj}), 1:num_rep, 'UniformOutput', false );
        is_empty_trial = cellfun(@(X) isscalar(X) && (0 == X), trials{ii,jj}, 'UniformOutput', true);
        trials{ii,jj} = trials{ii,jj}(~is_empty_trial);
    end
end



%% VERBOSE mode; write on command-line
if verbose
    fprintf('\n--------------------------------\n');
    fprintf('--> SESSION: %s\n', session_str);
    fprintf('--------------------------------\n');
    
    % Show the structure of the recorded data
    fprintf('--> S.t:\n');
    col_name_str = arrayfun(@(X) sprintf('Dist%d', X), S.outerSeq.master.values, 'UniformOutput', 0);
    row_name_str = arrayfun(@(X) sprintf('Revb%d', X), S.innerSeq.master.values, 'UniformOutput', 0);
    col_name_str = ['Revb', col_name_str];
    dummy_table = cell2table([row_name_str(:), S.t],...
        'VariableNames', col_name_str);
    disp(dummy_table)
    fprintf('\n');
    fprintf('\n--> syncchan:  %d\n', syncchan);
    fprintf('--> spikechan: %d\n', spikechan);
    %fprintf('--> All available spikechan(s): [%s]\n', num2str(all_spikechan, '%d '));
    %assert( logical(nnz(meas_st.spikechan==all_spikechan)),...
    %    '--> ERROR: the SPIKECHAN that you chose does NOT exist in this measurement!');

    %I_str = sprintf('%s: %d%%', S.innerSeq.master.var, 10*S.innerSeq.master.values(I));
    %J_str = sprintf('%s: %.1f m', S.outerSeq.master.var, 0.1*S.outerSeq.master.values(J));

    %fprintf('--> %s\n', I_str);
    %fprintf('--> %s\n', J_str);

end












