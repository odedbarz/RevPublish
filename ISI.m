function [H, bintimes, spikecount] = ISI(t, ch, t_start, t_end, varargin)
%
%   function [hdata, bintimes, spikecount, hc]  = ISI(t_start, t_end,
%       [t], [ch], [spikechan], [syncchan], [isi_order], [verbose]);
%
%
% * (First-order) interspike interval (ISI) histograms display the distribution of time intervals between consecutive spikes.  
% * They are most useful for spontaneous activity, and for periodic stimuli such as pure tones, where they reveal the presence of phase locking.
% * See also h.type == 'Interval' in hcompute.m at @histogram.
%

%% Parse the input
p = inputParser;

addRequired(p, 't', @iscell);
addRequired(p, 'ch', @iscell);
addRequired(p, 't_start', @isnumeric);
addRequired(p, 't_end', @isnumeric);


addOptional(p, 'binwidth', 1, @isnumeric);      % (ms)
addOptional(p, 'nbins', 200, @isnumeric);
addOptional(p, 'spikechan', 1, @isnumeric);
addOptional(p, 'syncchan', 0, @isnumeric);
addOptional(p, 'isi_order', 1, @isnumeric);
addOptional(p, 'verbose', 0, @isnumeric);
addOptional(p, 'fignum', 0, @isnumeric);

parse(p, t, ch, t_start, t_end, varargin{:});

nbins     = p.Results.nbins;
spikechan = p.Results.spikechan;
syncchan  = p.Results.syncchan;
binwidth  = p.Results.binwidth;
isi_order = p.Results.isi_order;
verbose   = p.Results.verbose;
fignum    = p.Results.fignum;


%%
assert(0 < t_end - t_start, '--> Error in [ana_c->ISI.m]: 0 >=  t_end - t_start!!');

h = set( histogram,...
    'Type', 'Interval',...
    'Size', [1, nbins],...      % [1, histSize], ...
    'BinWidths', binwidth, ...
    'Offsets', 0, ...           % [ms] when to start the ISI; !! note: this isn't the Gate!
    'Order', isi_order,...
    'SyncChan', syncchan, ...
    'SpikeChan', spikechan ...
);

% Set time window (a Gate) for the ISI histogram: 
g.type = 'PST';        % {'time', 'PST'}
g.delay = t_start;
g.width = t_end - t_start;
h = set(h, 'GateType', g.type, 'GateDelay', g.delay, 'GateWidth', g.width);

% get the desired results:
hobj       = cellfun(@(T,CH) hcompute(h,T,CH), t, ch, 'UniformOutput', false);         % histograms's handles
hdata      = cellfun(@(HIST) double(HIST)', hobj, 'UniformOutput', false);                 % the PSTs cell
bintimes   = cellfun(@(HIST) get(HIST, 'BinTimes'), hobj, 'UniformOutput', false);
spikecount = cellfun(@(HIST) get(HIST,'SpikeCount'), hobj, 'UniformOutput', true);        

% Arrange the data
%hdata      = cellfun(@(x) x(:), hdata, 'UniformOutput', false );
bintimes   = bintimes{1}(:);
H          = [hdata{:}];    % cells into one matrix

% Remove empty columns
H = H(:, 0<spikecount);


%% DEBUG MODE
if verbose
    fprintf('\n=== ISI Histogram: ===\n');
    fprintf('Type:       %s\n', 'Interval');
    fprintf('Size:       %g\n', nbins);
    fprintf('BinWidths:  %g ms\n', binwidth);
    fprintf('Hist Time:  %g ms\n', nbins * binwidth);
    fprintf('Offsets:    %g\n', t_start);    
    fprintf('ISI order:  %g\n', isi_order);
    fprintf('spikechan:  %g\n', spikechan);
    fprintf('syncchan:   %g\n', syncchan);    
    fprintf('\n');  
    fprintf('=== Gate: ===\n');
    fprintf('GateType:   %s\n', g.type);
    fprintf('GateDelay:  %g\n', g.delay);
    fprintf('GateWidth:  %g\n', g.width);
    fprintf('\n')
    
    fprintf('spikecount'); 
    disp( spikecount(:)' );
    fprintf('\n');  
end

if ~isempty(fignum) && fignum > 0
    figure(fignum);         
    clf;
    %subplot(2,1,1);
    plot(spikecount(:), '.--');
    xlabel('trial number');
    ylabel('Spikes');
    title('Spike Count');

    figure(1+fignum);         
    clf;    
    y_bias = max(H(:));
    for ii = 1:size(H,2)
        plot(bintimes, H(:,ii)/y_bias + ii);
        hold on       
    end
    hold off
    aux.hline(1:size(H,2), 'LineWidth', 0.5);
    set(gca, 'YTickLabel', 1:size(H,2))
    
    xlabel('Time [ms]');
    ylabel('trial number');
    title('ISI');
end


%% DEBUG
% Test the ISI histogram
%{
    N = 2e3;
    spikechan = 1;
    period_1 = 100;
    t = [0, find([0 == mod(1:N/2,period_1)])];
    period_2 = 70;
    t = [t, N/2 + find([0 == mod(1:N/2,period_2)])];

    noise_sd = 4;
    t = t + randi(noise_sd,size(t)) - noise_sd/2;
    t = {t(:)};

    ch = spikechan * ~[1 == 1:length(t{1})];
    ch = {ch(:)};

    [ch{1}, t{1}]

    ISI(0, 2e3, t, ch, 'binwidth', 1, 'spikechan', spikechan);
%}
        








