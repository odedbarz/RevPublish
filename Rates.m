function [rates, sem, hdata, hc] = Rates(t_start, t_end, t, ch, varargin)
%
%   [rates, sem, hdata, hc] = Rates(t_start, t_end, t, ch, [spikechan], [syncchan])
%
% Output:
%   * pst (=spikes per trial)
% 
%
% Notes:
% * PST Raster (a.k.a. dot raster) histograms show the distribution of spike times relative to stimulus onset along the horizontal axis, 
%   and for each stimulus presentation along the vertical axis.
% * A PST histogram is the column-by-column sum of the PST raster.
% * PST rasters give a nice visual display of the raw spike data very much as they are acquired during an experiment, 
%   and are useful for assessing the stability of the neural response over many repetitions of the same stimulus.
% * Appropriate bin widths for PST rasters are the same as for PST histograms.
%
% Oded Barzelay, 2/20/2018
%

%% Parse the input
p = inputParser;

addRequired(p, 't_start', @isnumeric);
addRequired(p, 't_end', @isnumeric);
addRequired(p, 't', @iscell);
addRequired(p, 'ch', @iscell);

addOptional(p, 'spikechan', 1, @isnumeric);
addOptional(p, 'syncchan', 0, @isnumeric);

parse(p, t_start, t_end, t, ch, varargin{:});

t_start   = p.Results.t_start;
t_end     = p.Results.t_end;
t         = p.Results.t;
ch        = p.Results.ch;
spikechan = p.Results.spikechan;
syncchan  = p.Results.syncchan;

%%
avg_dim = 1;    % the dimension to average over (usually, dim==1 is the F0)
% rates_dim = 2;  % the dimension along which to calc the rates (usually, dim==2 is the dB SPL)

% number of trials in each run:
ntrials = cellfun( @(CH) sum(syncchan == CH), ch, 'UniformOutput', false );

% Remove empty cells
n_trials= [ntrials{:}];
ntrials = ntrials(n_trials>0);
t       = t(n_trials>0);
ch      = ch(n_trials>0);

% Create a histograms object
h1 = set( histogram,...
    'Type', 'PST Raster',...
    'BinWidths', [1 t_end-t_start],...  % NOTE: the TOTAL width of the histogram is t_end = (t_end - t_start) + t_start
    'Offsets', [0 t_start],...
    'SyncChan', syncchan,...
    'SpikeChan', spikechan...
);


% Create a cell of histograms
hc = cellfun(@(X) h1, cell(size(t)), 'UniformOutput', false);  

% Set the Size for each of the histograms
hc = cellfun(@(H,NTRIALS) set(H,'Size',[NTRIALS, 1]), hc, ntrials, 'UniformOutput', false); 

% compute:
hc    = cellfun(@(H,T,CH) hcompute(H,T,CH), hc, t, ch, 'UniformOutput', false);   	% histograms's handles
hdata = cellfun(@(HIST) double(HIST), hc, 'UniformOutput', false);                  % the PSTs cell


%% Clac Rates & SEM
t_stim_eff_s = units.ms2sec( t_end - t_start );         % effective stimulus duration

% RATES (spikes/sec)
rates = cellfun( @(HDATA) mean(HDATA, avg_dim)/t_stim_eff_s, hdata, 'UniformOutput', true ); 

% SEM (standard error of the mean)
sem = cellfun( @(HDATA,NTR) std(HDATA, avg_dim)/( sqrt(NTR)*t_stim_eff_s ), hdata, ntrials, 'UniformOutput', true ); 























