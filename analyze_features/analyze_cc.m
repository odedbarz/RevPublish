%
% analyze_cc.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions.
%
%

clc
fignum = 10;
verbose = 1;

setup_environment('../');


analyze_setup;


%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 5;





%% Plot ALL CCt & CCt2 
figh = figure(fignum);
clf;

sp     = 1;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
idx_sp = idx_fun(sp);      % indices; time indices for speaker SP
drr_k  = 5;      % 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}

% Get the speaker's sex
if contains(tbl_metadata.fn(sp), '_M')
    sp_sex = 'Male';
else
    sp_sex = 'Female';
end

% Extract chunks\speakers
CCtk  = CCt(idx_sp, drr_k);
CCtk2 = CCt2(idx_sp, drr_k);

Sdry = spec_st.Sft{drr.ordered(1)};
Sdrr = spec_st.Sft{drr.ordered(drr_k)};
Sest = data.Sest;


%
tidx = t(idx_sp);       % (sec)
tidx = tidx-tidx(1);    % (sec)
plot(tidx, [CCtk, CCtk2], 'LineWidth', linewidth);
set(gca, 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
title(sprintf('%s Speaker (%d %ss)', sp_sex, n_units, data_type), 'FontSize', fontsize_big);

axis tight

legend( {sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}),...
    sprintf('CC$_t$($S_{%s}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}, drr_labels{drr_k}) },...
    'Location', 'southeast',...
    'FontSize', fontsize_big ...
);


% Set same positions for all figures
% pos_fig = [143, 330, 1610, 453];    % HP @ HOME
pos_fig = [143, 103, 1727, 680];    % HP @ MEE
set(figh, 'Position', pos_fig);






%% Plot a SAMPLE of the CCt & CCt2 with 3 highlight regions
%
figh = figure(2+fignum);
clf;

tidx = t(idx_sp);    % (sec)
tidx = tidx-tidx(1);    % (sec)
plot(tidx, [CCtk, CCtk2], 'LineWidth', linewidth);
set(gca, 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
title(sprintf('%s Speaker (%d %ss)', sp_sex, n_units, data_type), 'FontSize', fontsize_big);
axis tight

% legend('CC$_t$(Dry-Reconst.)', sprintf('CC$_t$(%s-Reconst.)', drr_labels{drr_k}),...
%     'Location', 'southeast', 'FontSize', fontsize_big);
legend( {sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}),...
    sprintf('CC$_t$($S_{%s}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}, drr_labels{drr_k}) },...
    'Location', 'southeast',...
    'FontSize', fontsize_big ...
);


%
face_alpha = 0.5;

% Time intervals fortidx the AREA
% '########## DEBUG ############'
% % ** sp = 10
% t_scope      = [1.31, 1.41];   % 1'st segment; $Dry > -8.2 dB$
% t_scope(2,:) = [1.42, 1.535];  % 2'nd segment; $Dry < -8.2 dB$
% t_scope(3,:) = [1.55, 1.65];   % 3'th segment; $Dry ~ %s$
% % ** sp = 1
% t_scope      = [0.546, 0.69];   % 1'st segment; $Dry > -8.2 dB$
% t_scope(2,:) = [0.705, 0.79];   % 2'nd segment; $Dry < -8.2 dB$
% t_scope(3,:) = [0.965, 1.1];    % 3'th segment; $Dry ~ %s$
% % ** sp = 9
% t_scope      = [0.535, 0.600];   % 1'st segment; $Dry > %s$
% t_scope(2,:) = [0.705, 0.79];   % 2'nd segment; $Dry < %s$
% t_scope(3,:) = [0.965, 1.1];    % 3'th segment; $Dry ~ %s$
% ** sp = 1
t_scope      = [2.84, 2.905];   % 1'st segment; $Dry > -8.2 dB$
t_scope(2,:) = [2.745, 2.815];  % 2'nd segment; $Dry < -8.2 dB$
t_scope(3,:) = [2.565, 2.635];   % 3'th segment; $Dry ~ %s$


hold on
area_h = area(t_scope(1,:), [1, 1]);
area_h(2) = area(t_scope(2,:), [1, 1]);
area_h(3) = area(t_scope(3,:), [1, 1]);
hold off

area_h(1).DisplayName = sprintf('$S_{dry} > S_{%s}$', drr_labels{drr_k});
area_h(1).FaceAlpha = face_alpha;
area_h(1).FaceColor = aux.rpalette('new03');
area_h(1).EdgeColor = 'k';

area_h(2).DisplayName = sprintf('$S_{dry} < S_{%s}$', drr_labels{drr_k});
area_h(2).FaceAlpha = face_alpha;
area_h(2).FaceColor = aux.rpalette('new04');
area_h(2).EdgeColor = 'k';

area_h(3).DisplayName = sprintf('$S_{dry} \\approx S_{%s}$', drr_labels{drr_k});
area_h(3).FaceAlpha = face_alpha;
area_h(3).FaceColor = aux.rpalette('new05');
area_h(3).EdgeColor = 'k';

% Set same positions for all figures
pos_fig = [113         228        1711         621];
set(figh, 'Position', pos_fig);


ylim(gca, [-2.0, 1.1]);
xlim(gca, [2.25, 2.98]);




%% IMAGESC of the previous 3 highlight regions to compare the Dry, Reconst., & DRR
%
clear figh ax
hz = 1e-3;     % {1.0, 1e-3}   % shows Hz or kHz
if (hz == 1e-3), hz_units='(kHz)'; else  hz_units='(Hz)'; end

Sdry_ = zca(Sdry(:,idx_sp));
Sest_ = zca(Sest(:,idx_sp, drr_k));
Sdrr_ = zca(Sdrr(:,idx_sp));

for k = 1:length(area_h)
    
    figh = figure(10+fignum + 2*k);
    clf;
    
    ax = subplot(1,3,1);
    img = imagesc(tidx, hz*f, Sdry_ );

    ax(2) = subplot(1,3,2);
    img(2) = imagesc(tidx, hz*f, Sest_ );

    ax(3) = subplot(1,3,3);
    img(3) = imagesc(tidx, hz*f, Sdrr_ );

    % Set the x-axis (time)
    linkaxes(ax);
    xlim(ax(1), t_scope(k,:));

    % set the color range for all imagesc
    min_cmap = 0;
    max_cmap = median( arrayfun(@(x) max(x.CData(:)), img) );
    arrayfun(@(AX) caxis(AX, [min_cmap, max_cmap]), ax );
    
    set(ax(2:3), 'YTickLabel', '');
    % set(ax(2:3), 'XTickLabel', '');
    set(ax, 'YDir', 'normal');
    set(ax, 'FontSize', fontsize);
    arrayfun(@(X) colormap(X, 'jet'), ax );

    ylabel(ax(1), ['Frequency ', hz_units], 'FontSize', fontsize_big);
    xlabel(ax(2), 'Time (sec)', 'FontSize', fontsize_big);

    title(ax(1), sprintf('$S_{dry}$ (%d %ss)', n_units, data_type), 'FontSize', fontsize_big);
    title(ax(2), sprintf('$\\hat{S}_{%s}$', drr_labels{drr_k}), 'FontSize', fontsize_big);
    title(ax(3), sprintf('$S_{%s}$', drr_labels{drr_k}), 'FontSize', fontsize_big);

    % Set same positions for all figures
    set(figh, 'Position', pos_fig);

    % move all axes left
    for nn = 1:length(ax)
        pos = get(ax(nn), 'Position');
        pos_new = pos;
        pos_new(1) = pos_new(1) + 0.05;
        pos_new(2) = pos_new(2) + 0.075;
        pos_new(4) = pos_new(4) - 0.15;
        set(ax(nn), 'Position', pos_new);
    end

    % set(fig, 'Color', h1.FaceColor);
    label_h = aux.abc(ax(1), 'location', 'northwest outside', ...
        'outside_xbias', 0.75 ,'outside_ybias', 0.12);
    label_h.BackgroundColor = area_h(k).FaceColor;
    label_h.FaceAlpha = face_alpha;
    label_h.EdgeColor = 'k';
    label_h.Color = 'none';
    % label_h.String = '';
    
    drawnow;

end






%% Plot \Delta CC_t
% This plot shows that the mean CC are in favor of the -8.2 dB-vs.-Reconst than the 
% Dry-vs.-Reconst reconstruction.
%{
figure(25 + fignum);
clf;
n_bins = 50;
%  ('count', 'probability', 'countdensity', 'pdf', 'cumcount', or 'cdf')
area_h = histogram(CCtk(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', 'CC$_t$(Dry vs. Reconst.)');
hold on
area_h(2) = histogram(CCtk2(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', sprintf('CC$_t$(%s vs. Reconst.)', drr.labels{drr.ordered(drr_k)}));
hold off
set(gca, 'FontSize', fontsize);

% aux.vline(nanmean(CCtk(:)-CCtk2(:)), 'LineStyle', ':' , 'displayName', 'mean()');
% legend('Location', 'northwest', 'FontSize', fontsize_big);
title(sprintf('Speaker %d, CC$_t$(Dry vs. Reconst.) vs CC$_t$(DRR vs. Reconst.)', sp), 'FontSize', fontsize_big);
ylabel('Count', 'FontSize', fontsize_big);
xlabel('CC', 'FontSize', fontsize_big);

mx = max([area_h(1).Values, area_h(2).Values]);
m1 = nanmean(area_h(1).Data);
m2 = nanmean(area_h(2).Data);

hold on
plot(m1, mx, 'v',...     
    'MarkerFaceColor', aux.rpalette('new01'),...
    'MarkerEdgeColor', 'none');
plot(m2, mx, 'v',...
    'MarkerFaceColor', aux.rpalette('new02'),...
    'MarkerEdgeColor', 'none');
hold off
legend(area_h, 'Location', 'northwest');
%}









%% PLOT Spectrogram + F0
figh = figure(40 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
ax = subplot(2,1,1);
[ax, surf_h] = spec.plot_spectrogram(pydata.spec.t, hz*pydata.spec.f, pydata.spec.Sdb,...
    'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', true);
set(ax, 'XTickLabel', '');
hold on
plot(t, log2(hz*F0i), 'w', 'LineWidth', 5);
hold off
title('Stimulus and its F0s');
ylabel(ax, 'Frequency (kHz)', 'FontSize', fontsize_big);
axis tight

%
ax(2) = subplot(2,1,2);
plot(ax(2), t, [CCt(:,n_drr), CCt2(:,n_drr)], 'LineWidth', linewidth);
set(ax(2), 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
axis tight
% hold on
% plot(t, F0i/max(F0i), 'LineWidth', linewidth);
% hold off
legend( {sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}),...
    sprintf('CC$_t$($S_{%s}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}, drr_labels{drr_k}) },...
    'Location', 'southeast',...
    'FontSize', fontsize_big ...
);
linkaxes(ax, 'x');
ax(1).Position(3:4) = ax(2).Position(3:4);
xlim([3.0, 12.0]);





%%
% Voiced regions are more robust for reverberation.
figure(45 + fignum);
clf;

y1 = CCt( idx_F0i, drr_k);
y2 = CCt(~idx_F0i, drr_k);

normalization = 'probability';  % {'pdf', 'probability', 'count', 'countdensity'}
h = histogram(y1, 50,...
    'Normalization', normalization,...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Voiced', drr_labels{drr_k}) ...
);
v_height = 0.8*max([h.Values]);
hold on
h(2) = histogram(y2, 50,...
    'Normalization', normalization, ...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Unvoiced', drr_labels{drr_k}) ...
);
plt = plot( nanmean(y1), v_height, 'v', ...
    'MarkerSize', markersize,...
    'MarkerFaceColor', aux.rpalette('new01'),...
    'MarkerEdgeColor', 'none' );
plt(2) = plot( nanmean(y2), v_height, 'v', ...
    'MarkerSize', markersize,...    
    'MarkerFaceColor', aux.rpalette('new02'),...
    'MarkerEdgeColor', 'none' );
hold off
set(gca, 'FontSize', fontsize);
xlabel('CC$_t$', 'FontSize', fontsize_big);
ylabel('Pr', 'FontSize', fontsize_big);
legend(h, 'Location', 'northwest');







%% Plot \Delta CC_t
% This plot shows that the mean CC are in favor of the -8.2 dB-vs.-Reconst than the 
% Dry-vs.-Reconst reconstruction.
% %{
figure(60 + fignum);
clf;
n_bins = 50;
%  ('count', 'probability', 'countdensity', 'pdf', 'cumcount', or 'cdf')
area_h = histogram(CCt(:)-CCt2(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', 'CC$_t$(Dry vs. Reconst.) - CC$_t$(DRR vs. Reconst.)');
set(gca, 'FontSize', fontsize);

title(sprintf('Speaker %d, CC$_t$(Dry vs. Reconst.) - CC$_t$(DRR vs. Reconst.)', sp), 'FontSize', fontsize_big);
ylabel('Count', 'FontSize', fontsize_big);
xlabel('$\Delta CC_t$', 'FontSize', fontsize_big);

% Fi = 
%      General model Gauss1:
%      Fi(x) =  a1*exp(-((x-b1)/c1)^2)
x = area_h.BinEdges(1:end-1) + diff(area_h.BinEdges);
y = area_h.Values;
Fi = fit(x', y', 'gauss1');

m0 = nanmean(CCtk(:)-CCtk2(:));

aux.vline(Fi.b1, 'Color', 'r');

x_ = interp1(1:length(x), x, 1:0.1:length(x));
hold on
plot(x_, Fi.a1*exp(-((x_-Fi.b1)./Fi.c1).^2), 'r');
plot(m0, 0.75*Fi.a1, 'rv', 'Color', aux.rpalette('new01'), 'MarkerFaceColor', aux.rpalette('new01'));
hold off
%}



%%
figure(63 + fignum);
clf;


% Indices of GOOD reconstructions
idx_robust = false(n_time, 1);
idx_robust( CCt(:, drr_k) - Fi.c1 >= CCt2(:, drr_k) ) = true;

% Keep only "goog" reconstructions
thr_good = 0.6;
idx_good = (CCt(:, drr_k)>=thr_good) .* (CCt2(:, drr_k)>=thr_good);

% n_prob = n_time;
N = nnz(idx_good);

plot(t, [CCt(:, drr_k), CCt2(:, drr_k)]);  
hold on
plot(t, F0i/max(F0i), 'LineWidth', 10);  
area_h = area(t, idx_robust);
hold off
set(gca, 'FontSize', fontsize);
area_h.FaceColor = aux.rpalette('new05');
area_h.FaceAlpha = face_alpha;
area_h.EdgeColor = 'none';
area_h.EdgeAlpha = 0;
xlim(gca, [3.0, 9.0]);
ylim(gca, [0.0, 1.0]);
ylabel('CC$_t$', 'FontSize', fontsize_big);
xlabel('Time (sec)', 'FontSize', fontsize_big);
legend( {sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}),...
    sprintf('CC$_t$($S_{%s}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}, drr_labels{drr_k}),...
    'F0/$F0_{max}$', ...
    sprintf('$S_{dry}$ + std $\\ge$ $S_{%s}$', drr_labels{drr_k}) },...
    'Location', 'southeastoutside',...
    'FontSize', fontsize_big ...
);


figure(65 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
num_F0_robust           = nnz(  idx_F0i .*  idx_robust .* idx_good )/N;
num_F0_not_robust       = nnz(  idx_F0i .* ~idx_robust .* idx_good )/N;
num_not_F0_robust       = nnz( ~idx_F0i .*  idx_robust .* idx_good )/N;
num_not_F0_not_robust   = nnz( ~idx_F0i .* ~idx_robust .* idx_good )/N;
M = [num_F0_robust,     num_F0_not_robust; ...
     num_not_F0_robust, num_not_F0_not_robust];

h = heatmap({'Dry + std > -8.2 dB', 'Dry + std < -8.2 dB'}, {'Harmonics', 'Non-Harmonics'}, M);
h.FontSize = fontsize; 






%% CCt + PHONEMES
figh = figure(70 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)

linewidth = 3

% AX #1 -- CCt
idx_sub = 4;
ax = subplot(20,1,1+idx_sub:20);
plot(ax, tidx, CCt(idx_sp,:), 'LineWidth', linewidth);
set(ax, 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
title(sprintf('%s Speaker (%d %ss)', sp_sex, n_units, data_type), 'FontSize', fontsize_big);
axis tight
legend_h = legend( drr.labels{drr.ordered}, 'Location', 'southeast' );
legend_h.FontSize = fontsize_big;

% AX #2 -- phonemes
ax(2) = subplot(20,1,1:idx_sub);
text_h = plot_phonemes(sp, 'ax', ax(2));

linkaxes(ax, 'x');

% ZOOM-IN
xlim_0 = 0.8;
xlim_1 = 2.2;
xlim([xlim_0, xlim_1]);
% Set same positions for all figures
% set(figh, 'Position', pos_fig);


 
 
%% SPECTROGRAM + F0 + PHONEMES
figh = figure(80 + fignum);
clf;
hz = 1e-3;

ti = nonzeros( idx_sp*(1e-3*binwidth) );

% AX #1 -- spectrogram
idx_sub = 4;
subplot(20,1,1+idx_sub:20);
[ax, surf_h] = spec.plot_spectrogram(ti-ti(1), hz*pyspec.f, pyspec.Sdbi(:, idx_sp),...
    'fontsize', fontsize, 'precision', 2, 'fignum', figh);

% Add F0s
hold on
plot(ax, ti-ti(1), log2(hz*F0i(idx_sp)), 'k', 'LineWidth', 5);
hold off


% AX #2 -- phonemes
ax(2) = subplot(20,1,1:idx_sub);
text_h = plot_phonemes(sp, 'ax', ax(2));

% Add the mean amplitudeto the phoneme's graph
mean_Si = mean(pyspec.Si(:,idx_sp));
hold on
plot(ax(2), ti-ti(1), mean_Si*(0.35/max(mean_Si)),...
    'Color', aux.rpalette('new01'), 'LineWidth', 5);
hold off

% ZOOM-IN
linkaxes(ax, 'x');
xlim_0 = 0.8;
xlim_1 = 2.2;
xlim([xlim_0, xlim_1]);


drawnow;
pos1 = get(ax(1), 'Position');
ax(2).Position(3) = pos1(3);


 
 
 










