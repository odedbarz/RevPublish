% plot_compare_scores_with_bars.m
%
% Compare scores of correlation coefficients (CCs) over various DRRs. The
% CCs are of the stimuli spectrograms vs. reconstructed spectrograms.
%
% This script is part of the analyze_reconst_scores.m script.

clc
% close all
% clear all

fignum = 51;
verbose = 1;

% Add the path of the StimViewerGUI GUI
% plot_compare_scores_with_bars_SAME_TRAINING_CONDITION.m
%

%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;


%% Load analyzed data
%     'scores2', 
%     'spec_st', 
%     'H', 
%     'splits', 
%     'n_neurons', 
%     'CC_broadband_envelope',...
%     'tbl_impale', 
%     'tbl_metadata'

% MUA
% %{
data_type = 'MUA';
analyze.path = '../.data/Analyze/';
analyze.filename = 'Analysis_MUA_(12-Oct-2020)_units(100)_bw(5)_fbands(30)_trainDRR(1-5).mat';
warning off
load([analyze.path, analyze.filename]);
warning on
%}




%% 

drr = get_DRR_list_and_indices;



%%
figure(0+fignum);
clf;

fontsize = 32;

aux.cprintf('String', 'Comparing Reconstructed Spectrograms of %s\n', data_type);

x = 1:drr.n_drr;
M = [scores.mu.CC(drr.ordered,end), scores_same_drr.mu.CC(drr.ordered,end)];
errbar = [scores.SE.CC(drr.ordered,end), scores_same_drr.SE.CC(drr.ordered,end)];
h = bar(M); 
set(gca, 'XTick', 1:drr.n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));

dx = 0.2*h(1).BarWidth;
hold on
er = errorbar([(1:5)'-dx, (1:5)'+dx], M, errbar);
hold off
for k = 1:length(er)
    er(k).Color = [0 0 0];                            
    er(k).LineStyle = 'none'; 
    er(k).LineWidth = 2;
end
ylim([0.0, 1.0]);
legend(h, {'To dry', 'To DRR'}, 'Location', 'northeast');
ylabel('CC($S(f,t)$-to-$\hat{S}(f,t)$)');
xlabel('DRR');
% aux.ctitle( sprintf('Comparing Reconstructed Spectrograms of %s', finfo.type),...
%     'to Dry and Other DRR Conditioned Stimuli');
title(sprintf('%d %ss', n_neurons, data_type ));
set(gca, 'FontSize', fontsize);

% #### set the FIGURE position ####
set(gcf, 'Position', [94, 2, 1192, 994])


