% plot_example_units_vs_time.m
%
% Plot few MUA units to show in the slide.
clc
fignum = 110;
verbose = 1;
addpath('../');
FigSetup;



%%
fn.load.path    = '../.data';
% fn.load.file    = 'data_MUA_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone).mat'; 
% fn.load.file = 'data_MUA_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)';
% fn.load.file = 'data_MUA_(23-Jul-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone)';
fn.load.file = 'data_SU_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone)';
fn.load.fullfile= fullfile( fn.load.path, fn.load.file );
data            = load(fn.load.fullfile);
spec_st         = data.spec_st;
duration_sec    = 1e-3*spec_st.duration_ms;  % (sec) stimulus duration to use
        
% Load the stimuli
dummy   = load(fullfile('../.data/stimulus', 'data_stimuli_duration(36_and_40)sec.mat'));
stim_st = dummy.stim_list{1};
assert(spec_st.duration_ms == stim_st.duration_ms);


drr = get_DRR_list_and_indices;

%%
figure(fignum);
clf;

fontsize = 24;
linewidth = 3;

% # of units to plot
N = 4;
% units = randperm(size(data.H,3), N);
units = [14, 16, 11, 2];    % SU selected 
% units = [23, 35, 104, 143];    % MUA selected

bw = data.spec_st.binwidth;
dt = 1e-3*bw;
t0 = 3.1;
n0 = max(1,ceil(t0/dt));
t1 = 6.25;
n1 = floor(t1/dt);
ns = n0:n1;
t = dt*(ns - ns(1));
nt = n1 - n0 +1;

Y = data.H(ns, drr.dry, units);

% Normalize for rates
Y = Y/(1e-3*bw);

ax = nan(1, N);
for k = 1:N
    ax(k) = subplot(N, 1, k);
    plot(t, Y(:,k), 'LineWidth', linewidth, 'Color', aux.rpalette(1));
    
    %xlim(ax(k), [t0, t1]-t0);        
    %ylim(ax(k), [0, 150]);
    
end

% axis(ax,'tight');
set(ax(1:end-1), 'XTickLabel', '', 'YTickLabel', '');
xlabel(ax(end), 'Time (sec)');
ylabel(ax(end),aux.ctitle('Firing Rate', '(spikes/sec)'));
set(ax, 'FontSize', fontsize);

