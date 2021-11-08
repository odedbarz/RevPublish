function rms_percentile_db = rms_running_window(y, Fs, varargin) 
%
%   function [rms_percentile_db, rms_percentile] = rms_running_window(y, Fs, 
%     percentile, win_size_ms, Nhist, fignum)
%
% Description:
% Finds the RMS value of a signal.
%

%% Parse the input
p = inputParser;

addRequired(p, 'y', @(x) isvector(x) && ~iscell(x));
addRequired(p, 'Fs', @isscalar);

addOptional(p, 'percentile', 95, @isscalar);
addOptional(p, 'win_size_ms', 40, @isscalar);
addOptional(p, 'Nhist', 100);
addOptional(p, 'overlap', 75, @(x) isscalar(x) && (x>0) && (x<=100));
addOptional(p, 'fignum', []);

parse(p, y, Fs, varargin{:});

percentile  = p.Results.percentile;            
win_size_ms = p.Results.win_size_ms;    
Nhist       = p.Results.Nhist;   
overlap     = p.Results.overlap;   
fignum      = p.Results.fignum;            


%%
% The size of the sliding window, in samples
win_size_smp = fix(1e-3*win_size_ms * Fs);

% The size of the overlap area
overlap_len = fix(overlap/100 * win_size_smp);
overlap_jump = max(win_size_smp - overlap_len, 1);

% Number of segments that the window slices over (discard the last samples
% at the end)
N = floor( (length(y) - win_size_smp)/overlap_jump) + 1;

% A vector that will hold all the RMS calculated by the sliding window
rms_v = zeros(1,N);

for kk = 1:N
    idx_k = (1:win_size_smp) + (kk-1)*overlap_jump;
    rms_v(kk) = rms( y(idx_k) );
end

rms_db_v = db(rms_v);   % 20*log10(*)

[counts, centers] = hist(rms_db_v, Nhist);
counts_sum        = cumsum(counts)/sum(counts);   % counts_sum: from 0 to 1
ind_percentile    = find( counts_sum >= percentile/100, 1, 'first' );
rms_percentile_db = centers(ind_percentile);



%% Plot results
if ~isempty(fignum)
    figure(fignum);
    clf;
    subplot(1,2,1);
    bar(centers, counts);
    hold on
    ywin = ylim;
    plot(rms_percentile_db*[1 1], [ywin(1), ywin(end)], '--k');
    hold off    
    title( sprintf('RMS (win size: %g ms)', win_size_ms) );
    xlabel('Root-mean-square levels');
    ylabel('Count');
    lgdh = legend('Count', sprintf('RMS$_{%g \\%%}$: %.1f dB',percentile, rms_percentile_db));
    set(lgdh, 'interpreter', 'latex');
    axis tight
    %
    subplot(1,2,2);
    bar(centers, counts_sum);
    grid on
    hold on
    plot([centers(1), centers(end)], percentile/100*[1, 1], '--k');
    ywin = ylim;
    plot(rms_percentile_db*[1 1], [ywin(1), ywin(end)], '--k');
    hold off
    title('Cumulative Sum of RMSs');
    lgdh = legend('$\Sigma rms_j$', sprintf('%g percentile', percentile/100));
    set(lgdh, 'interpreter', 'latex');
    xlabel('Root-mean-square levels');
    ylabel('Count (Normalized)');
    axis tight
end




















