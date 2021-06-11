function [h, plth, er] = plot_CC_S_vs_Sest(M, SE, CCs, title_str)
%
%   function [h, plth] = plot_CC_S_vs_Sest(M, ax)
%
% Ad hoc function for plotting Sdry/Sdrr vs. Sest bars.

import superbar.*

linewidth = 3;
markersize = 32;
fontsize = 24;
fontsize_big = 32;
fontsize_bigger = 42;
face_alpha = 1.0;

drr = get_DRR_list_and_indices;

% *** Sdry vs. Sest
h = bar(M); 
% h = superbar(M);
set(gca, 'FontSize', fontsize);
ylim([0.0, 1.2]);

set(gca, 'XTick', 1:drr.n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));
h(1).FaceAlpha = face_alpha;
h(2).FaceAlpha = face_alpha;

% CCs; add CC of stimuli
hold on
plth = plot(1:drr.n_drr, CCs, 'ks:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k' );

% Add error bars
dx = 0.18*h(1).BarWidth;
er = errorbar([(1:5)'-dx, (1:5)'+dx], M, SE, 'CapSize', 24);
for k = 1:length(er)
    er(k).Color = [0 0 0]; %h(k).FaceColor;                            
    er(k).LineStyle = 'none'; 
    er(k).LineWidth = linewidth;
end
hold off

legend(h, {'$S_{dry}$ vs. $\hat{S}_{drr}$', '$S_{drr}$ vs. $\hat{S}_{drr}$'},...
    'FontSize', fontsize_bigger);

ylabel('CC', 'FontSize', fontsize_big);
xlabel('DRR', 'FontSize', fontsize_big);

title(title_str, 'FontSize', fontsize_bigger);

