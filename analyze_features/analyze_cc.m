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
fontsize = 32;
fontsize_big = 42;
fontsize_bigger = 64;

markersize = 24;



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


% Extract data: 
CCt     = data.CCt;         % compares dry vs. est
CCt2    = data.CCt2;        % compares drr vs. est
splits  = data.splits;
spec_st = data.spec_st;
stim_st = data.stim_st;
tbl_data= data.tbl_data;
% n_splits= splits.n_splits;
n_bands = spec_st.n_bands;
n_time  = length(splits.idx);
% n_lags: # of samples along the x-axis (time) for each patch to use for comparison
patch_width = 1;    
Sest    = data.est_st.Sest;


% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
t           = spec_st.t;            % (sec)

drr     = get_DRR_list_and_indices; 
n_drr   = drr.n_drr;
idx_drr = drr.sortby(1:n_drr);
drr_labels = drr.labels(drr.ordered);


%% What to plot
sp = 5;
idx_x = sp == splits.idx;   % samples
k = idx_drr(5);

% Extract chunks\speakers
At  = CCt(idx_x, k);
At2 = CCt2(idx_x, k);
Sdry = spec_st.Sft{drr.ordered(1)};
Sdrr = spec_st.Sft{drr.ordered(k)};



%%
figh = figure(fignum);
clf;
tidx = t(idx_x);    % (sec)
tidx = tidx-tidx(1);    % (sec)
plot(tidx, [At, At2]);
legend('CC$_t$(Dry-Est)', sprintf('CC$_t$(%s-Est)', drr_labels{k}));
xlabel('Time (sec)');
ylabel('CC$_t$');
title(sprintf('%d %ss', n_units, data_type));
axis tight
ylim([-0.5, 2.0]);
xlim([0.55, 1.1]);

face_alpha = 0.5;

% Time intervals for the AREA
t_scope = [0.825, 0.925];
t_scope(2,:) = [0.74, 0.82];
t_scope(3,:) = [0.61, 0.725];

hold on
area_h = area(t_scope(1,:), [1, 1]);
area_h(2) = area(t_scope(2,:), [1, 1]);
area_h(3) = area(t_scope(3,:), [1, 1]);
hold off

area_h(1).DisplayName = sprintf('$Dry > %s$', drr_labels{k});
area_h(1).FaceAlpha = face_alpha;
area_h(1).FaceColor = aux.rpalette('new03');
area_h(1).EdgeColor = 'k';

area_h(2).DisplayName = sprintf('$Dry < %s$', drr_labels{k});
area_h(2).FaceAlpha = face_alpha;
area_h(2).FaceColor = aux.rpalette('new04');
area_h(2).EdgeColor = 'k';

area_h(3).DisplayName = sprintf('$Dry \\approx %s$', drr_labels{k});
area_h(3).FaceAlpha = face_alpha;
area_h(3).FaceColor = aux.rpalette('new05');
area_h(3).EdgeColor = 'k';

set(gca, 'FontSize', fontsize);

% Set same positions for all figures
% get(gcf, 'Position')
pos_fig = [113, 92, 1743, 655];
set(figh, 'Position', pos_fig);



%%
hz = 1e-3;     % {1.0, 1e-3}   % shows Hz or kHz
if (hz == 1e-3), hz_units='(kHz)'; else  hz_units='(Hz)'; end

for k = 1:length(area_h)
    
    figh = figure(10+fignum + 2*k);
    clf;
    
    ax = subplot(1,3,1);
    img = imagesc(tidx, hz*f, Sdry(:,idx_x) );

    ax(2) = subplot(1,3,2);
    img(2) = imagesc(tidx, hz*f, Sest(:,idx_x) );

    ax(3) = subplot(1,3,3);
    img(3) = imagesc(tidx, hz*f, Sdrr(:,idx_x) );

    linkaxes(ax);
    set(ax(2:3), 'YTickLabel', '');
    % set(ax(2:3), 'XTickLabel', '');
    set(ax, 'YDir', 'normal');
    set(ax, 'FontSize', fontsize);
    arrayfun(@(X) colormap(X, 'jet'), ax );

    ylabel(ax(1), ['Frequency ', hz_units]);
    xlabel(ax(2), 'Time (sec)');

    xlim(t_scope(k,:));

    title(ax(1), sprintf('Dry (%d %ss)', n_units, data_type));
    title(ax(2), 'Reconstructed');
    title(ax(3), 'DRR');

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
    % rectangle('Position', [20,20,100,100], 'FaceColor', [0 .5 .5])
    % annotation('rectangle',0.001*[20,20,100,100],'Color','red')
    label_h = aux.abc(ax(1), 'location', 'northwest outside', ...
        'outside_xbias', 0.75 ,'outside_ybias', 0.12);
    label_h.BackgroundColor = area_h(k).FaceColor;
    label_h.FaceAlpha = face_alpha;
    label_h.EdgeColor = 'k';
    label_h.Color = 'none';
    % label_h.String = '';

end

%%
% * If patch_width == 1 then there is no need to reshape columns...
% Rx = @(x) reshape(x, [n_bands, patch_width]);
% imagesc([Rx(At(:,nn)), Rx(At2(:,nn))]);
figure(20+fignum);
clf;
n_bins = 50;
area_h = histogram(At(:), n_bins);
hold on
h2 = histogram(At2(:), n_bins);
hold off















