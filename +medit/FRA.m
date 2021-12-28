function fra = FRA(S, spikechan, syncchan)
%
% function fra = FRA(S, spikechan, [syncchan])
%
% Input:
%   S         : (struct) an Impale's data structure of the FRA measurement.
%   spikechan : (1x1) the spike-channel requested number.
%   [syncchan]: (1x1) the syncchan, usually equals 0.
%
% Output:
%    fra: (struct) a structure that contains the CF, rates, and 
%     -fra.rates   : (NxM) the (dB x Frequency) rates.
%     -fra.CF      : (1x1) the characteristic frequency of the neuron.
%     -fra.graphics: (struct) holds data on the boundary line of neuron's most prominent region.

%

if 3 > nargin, syncchan = 0; end

gate        = find( arrayfun(@(ix) S.stimChans{ix}.IsActive, 1:length(S.stimChans)) );
assert(1 == length(gate), 'ERROR in [medit.FRA]: there are MORE THEN ONE one active gates!');

t_raisefall = S.stimChans{gate}.Gate.RiseFallTime;     % (ms)
t_delay     = S.stimChans{gate}.Gate.Delay;            % (ms) usually 0 for FRA
t_start     = t_delay + t_raisefall;                % (ms)
t_end       = S.stimChans{gate}.Gate.Width;            % (ms)



%% RATES: Set the spike times % keep the desired time intrval
[tps, nstim] = cellfun(@(T,CH) spiketools.pstimes(T,CH), S.t, S.ch, 'UniformOutput', false);
nstim = cell2mat(nstim);

idx.spikes = cellfun(@(CH) CH==spikechan, S.ch, 'UniformOutput', false);
idx.chans  = cellfun(@(CH) CH==syncchan, S.ch, 'UniformOutput', false);
idx.t      = cellfun(@(T) T>=t_start & T<=t_end , tps, 'UniformOutput', false);
idx.all    = cellfun(@(D1,D2,D3) (D1 & D2) | D3, idx.spikes, idx.t, idx.chans, 'UniformOutput', false);

% Get the selected spikes and sync channels
t_slc  = cellfun(@(T,IDX) T(IDX),    S.t, idx.all, 'UniformOutput', false);
ch_slc = cellfun(@(CH,IDX) CH(IDX), S.ch, idx.all, 'UniformOutput', false);

spike_count = cellfun(@(T,CH) numel(T(CH~=syncchan)), t_slc, ch_slc ,'UniformOutput', true);

% Total stimulus effective time
t_stim_eff_s = units.ms2sec( t_end - t_start );         % effective stimulus duration

% * mean fra_rates with respect to # of trials;
fra.rates = spike_count./nstim/t_stim_eff_s;



%% CF: calc the neuron's characteristic frequecny
[fra.CF, fra.graphics] = medit.calc_CF(S, fra.rates);



