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



%% Plot properties
fontsize = 24;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 24;
linewidth = 5;


%% Load the data
%              CC: [7200×5 double] 	% compares dry vs. est
%             CC2: [7200×5 double]	% compares drr vs. est
%             CCt: [7200×5 double]  % compares dry vs. est along the time domain
%            CCt2: [7200×5 double]  % compares drr vs. est along the time domain
%             PPt: [7200×5 double]
%            PPt2: [7200×5 double]
%     patch_width: 1
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [103×20 table]
data_type   = 'MUA';       % {'SU', MUA'}
switch data_type
    case 'SU'
        n_units = 103;
        
    case 'MUA'
        n_units = 241;   %[10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end



patch_width = 1;

fprintf('Loading the CC arrays...\n');
fn_path = '../_data/Analysis/';
fn_name = sprintf('analyzed_cc_%s_units(%d)_patchWidth(%d).mat', data_type, n_units, patch_width);
fn_fullfile = fullfile( fn_path, fn_name );
warning off
data = load(fn_fullfile);
warning on

% Print this to get more information about the data
if verbose
    aux.cprintf('Comments', 'Data Info:\n');
    aux.cprintf('Comments', data.info);
end

% Extract data: 
CC      = data.CC;         % compares dry vs. est
CC2     = data.CC2;        % compares drr vs. est
CCt     = data.CCt;         % compares dry vs. est
CCt2    = data.CCt2;        % compares drr vs. est
splits  = data.splits;
spec_st = data.spec_st;
stim_st = data.stim_st;
tbl_data= data.tbl_data;
n_bands = spec_st.n_bands;
n_time  = length(splits.idx);
% patch_width: # of samples along the x-axis (time) for each patch to use for comparison
patch_width = 1;    


% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
t           = spec_st.t;            % (sec)
drr         = get_DRR_list_and_indices; 
n_drr       = drr.n_drr;
drr_labels  = drr.labels(drr.ordered);





%% Plot ALL CCt & CCt2 
figh = figure(fignum);
clf;

sp = 1;
idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
drr_k = 5;      % 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}

% Extract chunks\speakers
CCtk  = CCt(idx_sp, drr_k);
CCtk2 = CCt2(idx_sp, drr_k);

Sdry = spec_st.Sft{drr.ordered(1)};
Sdrr = spec_st.Sft{drr.ordered(drr_k)};
Sest = data.Sest;


%
tidx = t(idx_sp);    % (sec)
tidx = tidx-tidx(1);    % (sec)
plot(tidx, [CCtk, CCtk2], 'LineWidth', linewidth);
set(gca, 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
title(sprintf('%d %ss', n_units, data_type), 'FontSize', fontsize_big);
axis tight

legend('CC$_t$(Dry-Est)', sprintf('CC$_t$(%s-Est)', drr_labels{drr_k}),...
    'Location', 'southeast', 'FontSize', fontsize_big);


% Set same positions for all figures
pos_fig = [143         330        1610         453];
set(figh, 'Position', pos_fig);






%% Plot a SAMPLE of the CCt & CCt2 with 3 highlight regions
% sp = 1;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
% drr_k = 5;      % 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}

% % Extract chunks\speakers
% CCtk  = CCt(idx_sp, drr_k);
% CCtk2 = CCt2(idx_sp, drr_k);

% Sdry = spec_st.Sft{drr.ordered(1)};
% Sdrr = spec_st.Sft{drr.ordered(drr_k)};
% Sest = data.Sest;


%
figh = figure(2+fignum);
clf;

tidx = t(idx_sp);    % (sec)
tidx = tidx-tidx(1);    % (sec)
plot(tidx, [CCtk, CCtk2], 'LineWidth', linewidth);
set(gca, 'FontSize', fontsize);
xlabel('Time (sec)', 'FontSize', fontsize_big);
ylabel('CC$_t$', 'FontSize', fontsize_big);
title(sprintf('%d %ss', n_units, data_type), 'FontSize', fontsize_big);
axis tight

legend('CC$_t$(Dry-Est)', sprintf('CC$_t$(%s-Est)', drr_labels{drr_k}),...
    'Location', 'southeast', 'FontSize', fontsize_big);


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

area_h(1).DisplayName = sprintf('$Dry > %s$', drr_labels{drr_k});
area_h(1).FaceAlpha = face_alpha;
area_h(1).FaceColor = aux.rpalette('new03');
area_h(1).EdgeColor = 'k';

area_h(2).DisplayName = sprintf('$Dry < %s$', drr_labels{drr_k});
area_h(2).FaceAlpha = face_alpha;
area_h(2).FaceColor = aux.rpalette('new04');
area_h(2).EdgeColor = 'k';

area_h(3).DisplayName = sprintf('$Dry \\approx %s$', drr_labels{drr_k});
area_h(3).FaceAlpha = face_alpha;
area_h(3).FaceColor = aux.rpalette('new05');
area_h(3).EdgeColor = 'k';

% Set same positions for all figures
pos_fig = [188 241 1479 608];
set(figh, 'Position', pos_fig);


ylim(gca, [-2.0, 1.1]);
xlim(gca, [2.25, 2.98]);




%% IMAGESC of the previous 3 highlight regions to compare the Dry, Est., & DRR
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

    title(ax(1), sprintf('Dry (%d %ss)', n_units, data_type), 'FontSize', fontsize_big);
    title(ax(2), 'Reconstructed', 'FontSize', fontsize_big);
    title(ax(3), 'DRR', 'FontSize', fontsize_big);

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
% This plot shows that the mean CC are in favor of the -8.2 dB-vs.-Est than the 
% Dry-vs.-Est reconstruction.
%
figure(20 + fignum);
clf;
n_bins = 50;
%  ('count', 'probability', 'countdensity', 'pdf', 'cumcount', or 'cdf')
area_h = histogram(CCt(:)-CCt2(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', 'CC$_t$(Dry-Est) - CC$_t$(DRR-Est)');
set(gca, 'FontSize', fontsize);

% aux.vline(nanmean(CCtk(:)-CCtk2(:)), 'LineStyle', ':' , 'displayName', 'mean()');
% legend('Location', 'northwest', 'FontSize', fontsize_big);
title(sprintf('Speaker %d, CC$_t$(Dry-Est) - CC$_t$(DRR-Est)', sp), 'FontSize', fontsize_big);
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
% plot(Fi.b1, 1.3*Fi.a1, 'rv', 'Color', 'r', 'MarkerFaceColor', 'r');
plot(m0, 0.75*Fi.a1, 'rv', 'Color', aux.rpalette('new01'), 'MarkerFaceColor', aux.rpalette('new01'));
hold off




%% Plot \Delta CC_t
% This plot shows that the mean CC are in favor of the -8.2 dB-vs.-Est than the 
% Dry-vs.-Est reconstruction.
%
figure(25 + fignum);
clf;
n_bins = 50;
%  ('count', 'probability', 'countdensity', 'pdf', 'cumcount', or 'cdf')
area_h = histogram(CCtk(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', 'CC$_t$(Dry vs. Est)');
hold on
area_h(2) = histogram(CCtk2(:), n_bins,...
    ...'Normalization', 'probability',...
    'displayName', sprintf('CC$_t$(%s vs. Est)', drr.labels{drr.ordered(drr_k)}));
hold off
set(gca, 'FontSize', fontsize);

% aux.vline(nanmean(CCtk(:)-CCtk2(:)), 'LineStyle', ':' , 'displayName', 'mean()');
% legend('Location', 'northwest', 'FontSize', fontsize_big);
title(sprintf('Speaker %d, CC$_t$(Dry-Est) vs CC$_t$(DRR-Est)', sp), 'FontSize', fontsize_big);
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





%% Loads the MAT file (LIBROSA)
% pydata = 
%   struct with fields:
% 
%     stim: [1×1 struct]
%        p: [1×1 struct]
%     harm: [1×1 struct]
%     spec: [1×1 struct]
mat_fn      = 'analyzed_librosa.mat';
mat_full_fn = fullfile( load.path_to_data('data'), 'Analysis', mat_fn);
pydata      = load( mat_full_fn );

% Extract parameters
spec = pydata.spec;
p   = pydata.p;

% sr = (pydata.stim.sr);
% duration_seconds = double(pydata.stim.duration_seconds);
% fmin = double(pydata.p.fmin);
% fmax = double(pydata.p.fmax);
% n_fft = (pydata.spec.n_fft);
% t = pydata.stim.t;
% y = double(pydata.stim.y);
% t_ = pydata.p.t_;
% F0 = pydata.p.F0;
% harm = pydata.harm;
% voiced_flag = pydata.p.voiced_flag;
% hop_length = double(pydata.p.hop_length);

% interpulate the F0 vector to the stimulus length
F0 = interp1( linspace(1,p.duration_seconds,p.nt_), p.F0, linspace(1,p.duration_seconds,n_time) )';


% PLOT Spectrogram + F0
figure(50 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)

imagesc(spec.t, hz*spec.f, spec.Sdb);
set(gca, 'YDir', 'normal');
colorbar;
ylim(hz*[0, 0.8e3]);
% set(gca, 'YScale', 'log')
set(gca, 'FontSize', fontsize);
hold on
plot(t, hz*F0, 'w');
hold off
title('Stimulus and its F0s');
ylabel('Frequency (kHz)', 'FontSize', fontsize_big);
xlabel('Time (sec)', 'FontSize', fontsize_big);



%%
% Indices of GOOD reconstructions
idx_robust = false(n_time, 1);
idx_robust( CCt(:, drr_k) + 2*Fi.c1/2 >= CCt2(:, drr_k) ) = true;

% Harmonic indices; get all indices of harmonics without consonants
idx_F0 = ~isnan(F0) .* F0>0;

% % Remove single elements (non-zeros surrounded by zeros)
% cond1 = idx_robust(1:end-2) .* idx_robust(3:end) .* ~idx_robust(2:end-1);     
% cond1 = [0; cond1; 0];
% plot([CCt(:, drr_k), CCt2(:, drr_k), idx_robust, cond1])


figure(52 + fignum);
clf;
% hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
% imagesc(spec.t, hz*spec.f, spec.Sdb);
% set(gca, 'YDir', 'normal');
% colorbar;
% hold on
% plot(t, hz*F0, 'w');
% plot(t, max(hz*F0)*idx_robust, 'w');
% hold off
plot([CCt(:, drr_k), CCt2(:, drr_k), F0/max(F0), idx_robust]);  '### DEBUG ###'


figure(55 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
num_F0_robust           = nnz(  idx_F0 .*  idx_robust )/n_time;
num_F0_not_robust       = nnz(  idx_F0 .* ~idx_robust )/n_time;
num_not_F0_robust       = nnz( ~idx_F0 .*  idx_robust )/n_time;
num_not_F0_not_robust   = nnz( ~idx_F0 .* ~idx_robust )/n_time;
M = [num_F0_robust,     num_F0_not_robust; ...
    num_not_F0_robust,  num_not_F0_not_robust];

h = heatmap({'Dry + std > -8.2 dB', 'Dry + std < -8.2 dB'}, {'Harmonics', 'Non-Harmonics'}, M);
h.FontSize = fontsize; 



nanmean(idx_F0 .* CCt(:, drr_k))
nanmean(idx_F0 .* CCt2(:, drr_k))
nanmean(~idx_F0 .* CCt(:, drr_k))
nanmean(~idx_F0 .* CCt2(:, drr_k))




%%
%{
figure(57 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
CCthr = 0.70;
num_F0_CCthr        = nnz(  idx_F0 .* CCt(:, drr_k)  > CCthr )/n_time;
num_not_F0_CCthr    = nnz( ~idx_F0 .* CCt(:, drr_k)  > CCthr )/n_time;
num_F0_CC2thr       = nnz(  idx_F0 .* CCt2(:, drr_k) > CCthr )/n_time;
num_not_F0_CC2thr   = nnz( ~idx_F0 .* CCt2(:, drr_k) > CCthr )/n_time;

bar([num_F0_CCthr, num_not_F0_CCthr; num_F0_CC2thr, num_not_F0_CC2thr]);



figure(55 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
num_F0_robust           = nnz(  idx_F0 .* CCt(:, drr_k)>CCthr   .*  idx_robust )/n_time;
num_F0_not_robust       = nnz( ~idx_F0 .* CCt(:, drr_k)>CCthr   .* ~idx_robust )/n_time;
num_not_F0_robust       = nnz(  idx_F0 .* CCt2(:, drr_k)>CCthr  .*  idx_robust )/n_time;
num_not_F0_not_robust   = nnz( ~idx_F0 .* CCt2(:, drr_k)>CCthr  .* ~idx_robust )/n_time;

M = [num_F0_robust,     num_F0_not_robust; ...
    num_not_F0_robust,  num_not_F0_not_robust];


bar(M);



%%
% Indices of GOOD reconstructions
CCthr = 0.60;
idx_robust = (CCt(:, drr_k) >= CCthr) .* (CCt2(:, drr_k) >= CCthr);


% Harmonic indices; get all indices of harmonics without consonants
idx_F0 = ~isnan(F0) .* F0>0;

% % Remove single elements (non-zeros surrounded by zeros)
% cond1 = idx_robust(1:end-2) .* idx_robust(3:end) .* ~idx_robust(2:end-1);     
% cond1 = [0; cond1; 0];
% plot([CCt(:, drr_k), CCt2(:, drr_k), idx_robust, cond1])


figure(52 + fignum);
clf;
% hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
% imagesc(spec.t, hz*spec.f, spec.Sdb);
% set(gca, 'YDir', 'normal');
% colorbar;
% hold on
% plot(t, hz*F0, 'w');
% plot(t, max(hz*F0)*idx_robust, 'w');
% hold off
plot([CCt(:, drr_k), CCt2(:, drr_k), F0/max(F0), idx_robust]);  '### DEBUG ###'


figure(55 + fignum);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
num_F0_robust           = nnz(  idx_F0 .*  idx_robust )/n_time;
num_F0_not_robust       = nnz(  idx_F0 .* ~idx_robust )/n_time;
num_not_F0_robust       = nnz( ~idx_F0 .*  idx_robust )/n_time;
num_not_F0_not_robust   = nnz( ~idx_F0 .* ~idx_robust )/n_time;
M = [num_F0_robust,     num_F0_not_robust; ...
    num_not_F0_robust,  num_not_F0_not_robust];

h = heatmap({'Dry & -8.2 dB > CCthr', 'Dry & -8.2 dB < CCthr'}, {'Harmonics', 'Non-Harmonics'}, M);
h.FontSize = fontsize; 


 
 

%%
% Indices of GOOD reconstructions
sf = spectral_flatness(spec.S);
sf = interp1( linspace(1,p.duration_seconds,p.nt_), sf, linspace(1,p.duration_seconds,n_time) )';

idx_robust = false(n_time, 1);
idx_robust( CCt(:, drr_k) + Fi.c1/2 >= CCt2(:, drr_k) ) = true;


histogram( (CCt(:, drr_k)-CCt2(:, drr_k) ) .* sf)
hold on
histogram(CCt2(:, drr_k) .* sf)
hold off
 
%}
 
 
 
 
 
 
 
 
 
 
 
 
 










