% plot_STRFs_vs_DRRs.m
%
%
clc;
fignum = 11;

addpath('../');
FigSetup;


%% Load the STRFs data of SUs
% whos
%   Name               Size                  Bytes  Class     Attributes
% 
%   H_valid         7200x6x103            35596800  double              
%   scores             1x1                  149976  struct              
%   slc                1x1                    1326  struct              
%   spec_st            1x1                15039301  struct              
%   splits             1x1                   58592  struct              
%   tbl_impale       437x20                 640191  table               
%   tbl_strf         103x5                  373423  table               
fn.path = '..\.data\STRFs\';
% fn.name = 'STRF_SU_(25-Aug-2020)_units(103)_bw(5)ms_algo(regression)_fbands(30)_lags(30)ms_cau(1)_trainDRR(3).mat';
fn.name = 'STRF_SU_(27-Aug-2020)_units(103)_bw(1)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(3)';
fn.file2load = fullfile(fn.path, fn.name);
clear data
data{1} = load(fn.file2load);

fn2.path = fn.path;
fn2.name = 'STRF_SU_(27-Aug-2020)_units(103)_bw(1)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(4)';
fn2.file2load = fullfile(fn2.path, fn2.name);
data{2} = load(fn2.file2load);
assert(data{2}.spec_st.binwidth == data{1}.spec_st.binwidth);

drr = get_DRR_list_and_indices;
%% Plot best frequency (BF) histogram
figure(fignum);
clf;

bf1         = 1e-3*data{1}.tbl_strf.bf;   % (kHz)
bf2         = 1e-3*data{2}.tbl_strf.bf;   % (kHz)
boundary    = 2.^(1.0/12) - 1;   % 1 OCTAVE
inv_bf      = ((1.0+boundary) >= bf1./bf2) & ((1.0-boundary) <= bf1./bf2);  % invariant BFs
units_inv_list = find(inv_bf)';
slc_idx_2_show = [4, 19, 17];
unit_idx    = units_inv_list(slc_idx_2_show);

fontsize = 24;
markersize = 34;


plot(bf1, bf2, '.', 'MarkerSize', markersize);
hold on
plot([min(xlim), max(xlim)], [min(ylim), max(ylim)], '--k');
nnz(inv_bf)
plth = plot(bf1(inv_bf), bf2(inv_bf), 'ok');            % plot the ONE-OCTAVE deviation
plth(2) = plot([0.2, 8], (1-boundary)*[0.2, 8], ':k');
plth(3) = plot([0.2, 8], (1+boundary)*[0.2, 8], ':k');
for q = 1:length(unit_idx)
    plth(3+q) = plot(bf1(unit_idx(q)), bf2(unit_idx(q)), 'o',...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', aux.rpalette(q+1),...
        'MarkerSize', 0.5*markersize);     % plot the selected unit
end
hold off
legend(plth(2), {'1 Octave Dev.'}, 'Location', 'northwest');
axis square
xlabel(aux.ctitle('Frequency (kHz)', 'Dry'));
ylabel(aux.ctitle('DRR(-8.2 dB)', 'Frequency (kHz)'));
set(gca, 'FontSize', fontsize);
title(sprintf('%d SUs', height(data{2}.tbl_strf)));
xlim(1e-3*[250, 8000]);
ylim(1e-3*[250, 8000]);



%% Plot the STRFs
figure(10+fignum);
clf;

f        = 1e-3*data{1}.spec_st.f;  % (kHz)
binwidth = data{1}.spec_st.binwidth;    % (ms)
fontsize = 24;

drr_ordered = [3, 4]

N = length(data);
M = length(unit_idx);
ax = nan(N, M);
for q = 1:M
    
    for k = 1:N
        %ax(k,q) = subplot(N,M,q+(k-1)*N);
        ax(k,q) = subplot(M,N,N*(q-1)+k);
        idx = drr_ordered(k);

        data{k}.obj.strf = data{k}.tbl_strf.strf{unit_idx(q)};
        bf_qk = 1e-3*data{k}.tbl_strf.bf(unit_idx(q));      % (kHz)
        data{k}.obj.plot_strf;
        colorbar('off');
        
        % Set the labels
        if M > q
            xlabel('');
        end
        if 1 == q
            title( sprintf('%s', drr.labels{idx}) );
        else
            title('');
        end
        
        if 1 < k
            ylabel('');
        else
            ylabel(aux.ctitle(sprintf('%.2f kHz', bf_qk), 'Frequency (kHz)'));
        end
        
    end
    
end

set(ax(:,1:end-1), 'XTickLabels', '');
% xlabel(ax(2:end,:), '');

set(ax, 'FontSize', fontsize);
set(ax, 'FontSize', fontsize);



%% Is the "invariant STRFs" keep their shape between DRRs?
n_units = height(data{1}.tbl_strf);
strf_diff = nan(n_units, 1);
for k = 1:n_units
    strf_diff(k) = norm( data{1}.tbl_strf.strf{k}- data{2}.tbl_strf.strf{k} );
end

% plot(strf_diff, '.');
% hold on
% plot(find(inv_bf), strf_diff(inv_bf), 'ok');
% hold off

h1 = histogram(strf_diff, 20, 'FaceColor', aux.rpalette('new01'));
hold on
h2 = histogram(strf_diff(inv_bf), 20, 'FaceColor', aux.rpalette('new02'));
hold off

