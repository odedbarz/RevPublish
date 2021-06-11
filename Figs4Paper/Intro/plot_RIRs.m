
clc
fignum = 10;
verbose = 1;

addpath('../../');
setup_environment;

analyze_setup;



%% Plot properties
fontsize = 58;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 6;


%%
fn = fullfile(load.path_to_data('Stimulus'), 'Spch_(36)sec/_TIMIT_meta.mat');
dummy = load(fn);
tbl_rir = dummy.Trir;


%%
figh = figure(fignum);
clf;

fs = tbl_rir.Fs(1);
x5 = tbl_rir.rir{end}(:,1);
% find(t>29, 1, 'first')
x5_early = nan*x5;
x5_early(202:2602) = x5(202:2602);    	% early reflections
x5_delta = nan*x5;
x5_delta(1:202) = x5(1:202);           % direct sound
x5(1:2602) = nan;

t = 1000*(0:length(y)-1)'/fs;

plth = plot(t, x5);
set(gca, 'FontSize', fontsize);
xlabel('Time (ms)');
ylabel('Amplitude');
xlim([-10, 300]);
ylim([-0.2, 1.1]);
% set(plth, 'LineWidth', 3);
title('RIR (-8.2 dB, Left Ear)');

hold on
plth(2) = plot(t, x5_early);
plth(3) = plot(t, x5_delta);
hold off
set(plth, 'LineWidth', linewidth);


