%
% plot_CCt.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions.
%
%

clc
fignum = 10;
verbose = 1;

setup_environment('../');


% Loads all analysis data 
analyze_setup;



%% Plot properties
fontsize = 36;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 6;



%% Analyze parameters
sp      = 10;             % speaker number to show

% Get the speaker's sex
if contains(tbl_metadata.fn(sp), '_M')
    sp_sex = 'Male';
else
    sp_sex = 'Female';
end

% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
idx_sp  = idx_fun(sp);	% indices; time indices for speaker SP
drr_k   = 5;             % 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}
tidx    = t(idx_sp);       % (sec)
tidx    = tidx-tidx(1);    % (sec)

Sdry = spec_st.Sft{drr.ordered(1)};
Sdrr = spec_st.Sft{drr.ordered(drr_k)};
% Sest = data.Sest;

% Energy as a function of time
ESdry = sum(Sdry);
ESdrr = sum(Sdrr);


% Get the speaker's sex
if contains(tbl_metadata.fn(sp), '_M')
    sp_sex = 'Male';
else
    sp_sex = 'Female';
end

% Extract chunks\speakers
CCtk  = CCt(idx_sp, drr_k);
CCtk2 = CCt2(idx_sp, drr_k);
CCtk3 = CCt3(idx_sp, drr_k);






%% SPECTROGRAM + F0 + PHONEMES + CCt, CCt2 & CCt3
figh = figure(fignum);
clf;
hz = 1e-3;

% CORRELATION
Cs = corrcoef([CCt(idx_nn,5), CCt2(idx_nn,5), CCt3(idx_nn,5)]);
fprintf('\nCs:\n');
disp(Cs);


ti = nonzeros( idx_sp*(1e-3*binwidth) );

% AX #1 -- Spectrogram
sub_0 = 9;
sub_1 = 11;
ax = subplot(20,1,2:sub_0);
% [~, surf_h] = spec.plot_spectrogram(ti-ti(1), hz*pyspec.f, pyspec.Sdbi(:, idx_sp),...
%     'fontsize', fontsize, 'precision', 2, 'fignum', figh);
[~, surf_h] = spec.plot_spectrogram(ti-ti(1), hz*f, spec_st.Sft{drr.dry}(:, idx_sp),...
    'fontsize', fontsize, 'precision', 2, 'fignum', figh);
set(ax(1), 'XTickLabel', '');
xlabel(ax(1), '');
ylabel(ax(1), 'Frequency (kHz)', 'FontSize', fontsize_big);
colorbar(ax(1));

% title_str = sprintf('%s Speaker Uttered "%s"', sp_sex, tbl_metadata.txt{sp}(9:end-1));
%title( title_str, 'FontSize', fontsize_big);
aux.ctitle(...
        sprintf('%s Speaker Uttered', sp_sex),...
        ['"', tbl_metadata.txt{sp}(9:end-2), '"'],...
    'FontSize', fontsize_big );



% AX #2 -- Phonemes
ax(2) = subplot(20,1,1+sub_0:sub_1);
text_h = plot_phonemes(sp, 'ax', ax(2));
set(text_h, 'FontSize', fontsize_big)
hold on
plth_ax2 = plot(tidx, 0.01*vc_nans(idx_sp), 'LineWidth', 20, 'Color', aux.rpalette('04'),...
    'DisplayName', 'Voiced Regions');
hold off
% set(ax(2), 'FontSize', fontsize_big);


% AX #3 -- CCts
ax(3) = subplot(20,1,sub_1+1:20);
%
plth = plot(tidx, CCtk,...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}));
hold on
plth(2) = plot(tidx, CCtk2,...
    'DisplayName', sprintf('CC$_t$($S_{%s}$ vs. $\\hat{S}_{%s}$)', drr_labels{drr_k}, drr_labels{drr_k}));
plth(3) = plot(tidx, CCtk3, 'k-',...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $S_{%s}$)', drr_labels{drr_k}) );
hold off
set(ax(3), 'FontSize', fontsize);
% set(plth(1:2), 'LineWidth', 2*linewidth);
% set(plth(3), 'LineWidth', linewidth);
set(plth, 'LineWidth', linewidth);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);

legend_h = legend(ax(3), 'Location', 'southeast' );
legend_h.FontSize = fontsize_big;
ylim([-1, 1]);



% ZOOM-IN
linkaxes(ax, 'x');
axis tight
% xlim_0 = 0.8;
% xlim_1 = 2.2;
% xlim([xlim_0, xlim_1]);

drawnow;
% pos1 = get(ax(1), 'Position');
% ax(1).Position(3) = pos1(3);
% ax(2).Position(3) = pos1(3);
% ax(3).Position(3) = pos1(3);
ax(1).Position([1,3]) = [0.09, 0.85];
ax(2).Position([1,3]) = [0.09, 0.85];
ax(3).Position([1,3]) = [0.09, 0.85];







%%
% Voiced regions are more robust for reverberation.
figure(10 + fignum);
clf;

n_bins = 50;

y_voiced   = CCt( idx_vc, drr_k);
y_unvoiced = CCt(~idx_vc, drr_k);

y_voiced3   = CCt3( idx_vc, drr_k);
y_unvoiced3 = CCt3(~idx_vc, drr_k);

% Histograms of Sdry-vs-\hat{S}drr
normalization = 'probability';  % {'pdf', 'probability', 'count', 'countdensity'}
h = histogram(y_voiced, n_bins,...
    'FaceColor', aux.rpalette('new04'),...
    'Normalization', normalization,...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Voiced', drr_labels{drr_k}) ...
);
v_height = 0.8*max([h.Values]);
hold on
h(2) = histogram(y_unvoiced, n_bins,...
    'FaceColor', aux.rpalette('new05'),...    
    'Normalization', normalization, ...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Unvoiced', drr_labels{drr_k}) ...
);


% % Histograms (lines) of Sdry-vs-Sdrr
% [N, edges] = histcounts(y_voiced3, n_bins);
% N_ = N/sum(N);
% edges_ = edges(1:end-1) + diff(edges);
% plt = plot(edges_, N_, 'Color', aux.rpalette('01'), 'LineWidth', linewidth,...
%     'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $S_{%s}$) Voiced', drr_labels{drr_k}));
% [N, edges] = histcounts(y_unvoiced3, n_bins);
% N_ = N/sum(N);
% edges_ = edges(1:end-1) + diff(edges);
% plt(2) = plot(edges_, N_, 'Color', aux.rpalette('02'), 'LineWidth', linewidth,...
%     'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $S_{%s}$) Unvoiced', drr_labels{drr_k}));


% Add MARKERS
pltm = plot( nanmean(y_voiced), v_height, 'v', ...
    'MarkerSize', markersize,...
    'MarkerFaceColor', h(1).FaceColor,...
    'MarkerEdgeColor', 'none' );
pltm(2) = plot( nanmean(y_unvoiced), v_height, 'v', ...
    'MarkerSize', markersize,...    
    'MarkerFaceColor', h(2).FaceColor,...
    'MarkerEdgeColor', 'none' );
% pltm(3) = plot( nanmean(y_voiced3), 1.2*v_height, 'v', ...
%     'MarkerSize', markersize,...
%     'LineWidth', linewidth,...
%     'Color', aux.rpalette('new01'),...
%     'MarkerFaceColor', 'none' );
% pltm(4) = plot( nanmean(y_unvoiced3), 1.2*v_height, 'v', ...
%     'MarkerSize', markersize,...  
%     'LineWidth', linewidth,...
%     'Color', aux.rpalette('new02'),...    
%     'MarkerFaceColor', 'none' );
hold off
set(gca, 'FontSize', fontsize_big);
xlabel('CC$_t$', 'FontSize', 1.2*fontsize_big);
ylabel('Pr', 'FontSize', 1.2*fontsize_big);
% legend([h, plt], 'Location', 'northwest', 'FontSize', fontsize_big);
legend(h, 'Location', 'northwest', 'FontSize', 1.2*fontsize_big);

ylim([0, 0.2]);





%% COMMAND-LINE
% CCt_to_CCt3_dry = corrcoef(CCt(idx_nn,1), CCt3(idx_nn,1));
% CCt_to_CCt3_dry = CCt_to_CCt3_dry(1,2);

CCt_1to3 = corrcoef(CCt(idx_nn,drr_k), CCt3(idx_nn,drr_k));
CCt_1to3 = CCt_1to3(1,2);

CCt_1to2 = corrcoef(CCt(idx_nn,drr_k), CCt2(idx_nn,drr_k));
CCt_1to2 = CCt_1to2(1,2);

% fprintf('CCt_1to3: %g\n', CCt_1to3);
fprintf('\n');
fprintf('CCt_1to3 : %.3f\n', CCt_1to3);
fprintf('CCt_1to2 : %.3f\n', CCt_1to2);






%% Dist. of ENERGY between VOICED & UNVOICED
% Voiced regions are more robust for reverberation.
%{
figure(40 + fignum);
clf;

n_bins = 50;

y_voiced   = CCt( idx_vc, drr_k);
y_unvoiced = CCt(~idx_vc, drr_k);

% y_voiced3   = ESdry( idx_vc)/max(ESdry);
% y_unvoiced3 = ESdry(~idx_vc)/max(ESdry);
y_voiced3   = ESdrr( idx_vc)/max(ESdrr);
y_unvoiced3 = ESdrr(~idx_vc)/max(ESdrr);

normalization = 'probability';  % {'pdf', 'probability', 'count', 'countdensity'}
h = histogram(y_voiced, n_bins,...
    'Normalization', normalization,...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Voiced', drr_labels{drr_k}) ...
);
v_height = 0.8*max([h.Values]);
hold on
h(2) = histogram(y_unvoiced, n_bins,...
    'Normalization', normalization, ...
    'DisplayName', sprintf('CC$_t$($S_{dry}$ vs. $\\hat{S}_{%s}$) Unvoiced', drr_labels{drr_k}) ...
);

[N, edges] = histcounts(y_voiced3, n_bins, 'Normalization', normalization);
N_ = N/sum(N);
edges_ = edges(1:end-1) + diff(edges);
plt = plot(edges_, N_, 'Color', aux.rpalette('01'), 'LineWidth', linewidth,...
    'DisplayName', sprintf('$E\\{S_{dry}\\}/max(E)$ Voiced', drr_labels{drr_k}));
[N, edges] = histcounts(y_unvoiced3, n_bins, 'Normalization', normalization);
N_ = N/sum(N);
edges_ = edges(1:end-1) + diff(edges);
plt(2) = plot(edges_, N_, 'Color', aux.rpalette('02'), 'LineWidth', linewidth,...
    'DisplayName', sprintf('$E\\{S_{dry}\\}/max(E)$ Unvoiced'));

title('Distribution of Energy between Voiced and Unvoiced');

% Add MARKERS
pltm = plot( nanmean(y_voiced), v_height, 'v', ...
    'MarkerSize', markersize,...
    'MarkerFaceColor', aux.rpalette('new01'),...
    'MarkerEdgeColor', 'none' );
pltm(2) = plot( nanmean(y_unvoiced), v_height, 'v', ...
    'MarkerSize', markersize,...    
    'MarkerFaceColor', aux.rpalette('new02'),...
    'MarkerEdgeColor', 'none' );
hold off
set(gca, 'FontSize', fontsize);
xlabel('CC$_t$', 'FontSize', fontsize_big);
ylabel('Pr', 'FontSize', fontsize_big);
legend([h, plt], 'Location', 'northwest', 'FontSize', fontsize_bigger);





%}














