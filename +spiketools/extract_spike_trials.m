function spks = extract_spike_trials(S, spikechan, I, J, syncchan)
%
%   function handles = set_tabs(filename)
%
% Input:
%   S: (struct) Impale's structure.
%   spikechan: (int) a spike-channel number to load.
%   I, J: indices to load from S.t & S.ch
%
% Output:
%   spks.spk_array: spikes cell array; contains just spikes, WITHOUT the
%              syncchan(s).
%   spks.tps: spike times relative to the onset.
%
%
% Description:
% Extracts the spike-times (tps; i.e., relative to the onset) of the desired 
% spikechan with their respective syncchan.
%

if 5 > nargin, syncchan = 0; end;

assert(isstruct(S) && isfield(S, 't') && isfield(S, 'ch'),...
    '--> [extract_spikes.m] <filename> MUST point into Impale''s structure!!');

assert( isfield(S.info,'reordered') && 1==S.info.reordered, ...
    '--> [extract_spikes.m]: The loaded Impale structure MUST be REORDERED (see the StimViewerGUI)!!!');

% Get a list of all spikechan 
all_spikechan = spiketools.get_all_spikechan(S);
assert(0 < nnz(all_spikechan == spikechan), '--> [extract_spikes.m] the requested SPIKECHAN isn''t available in S!!!');


%% Get the indices of the desired SPIKECHANs (plus their SPIKESYNC)
%[t, ch, ~] = extract_trials(S, I, J, spikechan, syncchan);
spk_idx = S.ch{I,J} == spikechan;   % get the spikechan spikes
trials = cumsum(0==S.ch{I,J});      % number all trials

% Select only the trial numbers of the desired SPIKECHAN (without the syncchan)
trials_of_spikechan = setdiff(unique(trials .* spk_idx), syncchan);

% Get all desired trials (spikechan + syncchan)
trials_idx = false(size(trials));
for kk = 1:length(trials_of_spikechan)
    trials_idx(trials == trials_of_spikechan(kk)) = true;
end

t  = S.t{I,J}(trials_idx);
ch = S.ch{I,J}(trials_idx);


%% Get spike-times in msec
[tps, nstim] = spiketools.pstimes(t, ch, syncchan);
if 0 == nstim
    spks = cell(1, nstim);
    return;
end

% just_spikes indices, without syncchan(s)
just_spikes = ch ~= syncchan;

trials = cumsum(ch == 0);
spk_array = arrayfun(@(I) tps(trials==I & just_spikes), 1:nstim, 'UniformOutput', false);


%% Save for output argument
spks.spk_array  = spk_array;
spks.tps        = tps;
spks.nstim      = nstim;
spks.trials_idx = trials_idx;
spks.t          = t;
spks.ch         = ch;




















