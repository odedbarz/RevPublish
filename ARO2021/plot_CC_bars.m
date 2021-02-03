% plot_CC_bars.m
%
% Plots error bars for Dry-vs-DDR reconstructions.


clc
fignum = 10;
verbose = 1;

setup_environment('../');


%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(3*fontsize); %64;

markersize = 36;
linewidth = 5;


%% Analyze
analyze_setup;

CCmu = squeeze( median(CC, 2) );
CCse = squeeze( mad(CC,[],2)/sqrt(size(CC,2)) );

CCmu2 = squeeze( median(CC2, 2) );
CCse2 = squeeze( mad(CC2,[],2)/sqrt(size(CC2,2)) );



%%
figure( 0 + fignum );
clf;

x = 1:n_drr;
M = [CCmu, CCmu2];
errbar = [CCse, CCse2];
h = bar(M); 
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));

dx = 0.2*h(1).BarWidth;
hold on
% ERROR BARs
er = errorbar([(1:5)'-dx, (1:5)'+dx], M, errbar,'CapSize',24);

% CCs
plth = plot(1:n_drr, CCs, 'ks:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k' );
hold off
set(gca, 'FontSize', fontsize_big);

for k = 1:length(er)
    er(k).Color = [0 0 0];                            
    er(k).LineStyle = 'none'; 
    er(k).LineWidth = 4;
end
ylim([0.0, 1.0]);
legend(h, {'$S_{dry}$ vs. $\hat{S}_{drr}$', '$S_{drr}$ vs. $\hat{S}_{drr}$'}, 'Location', 'northeastoutside', 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_bigger);
xlabel('DRR', 'FontSize', fontsize_bigger);
% aux.ctitle( sprintf('Comparing Reconstructed Spectrograms of %s', data_type),...
%     'to Dry and Other DRR Conditioned Stimuli');
title(sprintf('%d %s',n_units,data_type), 'FontSize', fontsize_bigger);









