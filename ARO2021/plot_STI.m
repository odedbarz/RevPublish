% plot_STI.m
%
% Plots STI scores

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
fn.path = '..\.data\STI\';

% Load the SU STI data
fn.name = 'STI_MUA_(03-Sep-2020)_units([100 150])_bw(5)_fbands(30).mat';
fn.file2load = fullfile(fn.path, fn.name);
mua = load(fn.file2load);

% Load the MUA STI data
fn.name = 'STI_SU_(03-Sep-2020)_units(103)_bw(5)_fbands(30).mat';
fn.file2load = fullfile(fn.path, fn.name);
su = load(fn.file2load);



%% Plot Speech Transmission Index (STI)
figure(fignum);
clf;

fontsize = 48;
markersize = 20;
linewidth = 3;

% The data matrix for the bar plot
M      = mean( su.sti.est(drr.ordered,:), 2 );
M(:,2) = mean( mua.sti.est(drr.ordered,:,1), 2 );   % x100 units
M(:,3) = mean( mua.sti.est(drr.ordered,:,2), 2 );   % x150 units

% The "test" probe (stimulus)
Mt = mean( su.sti.test(drr.ordered, :), 2);

SD      = std( su.sti.est(drr.ordered,:), [], 2 );
SD(:,2) = std( mua.sti.est(drr.ordered,:,1), [], 2 );   % x100 units
SD(:,3) = std( mua.sti.est(drr.ordered,:,2), [], 2 );   % x150 units

bar(M);
hold on
plot(Mt, 'sk', 'MarkerFaceColor', 'k',...
    'MarkerSize', markersize,...
    'LineWidth', linewidth);
% Add error bars for all 3 cases
for k = 1:3
    errorbar(-0.24+0.24*(k-1)+(1:5), M(:,k), SD(:,k), '.',...
        'LineWidth', linewidth, ...
        'Color', 'k'); %aux.rpalette(k));
end
hold off
ax = gca;
set(ax, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.ordered));
xlabel('DRR');
ylabel('STI');
set(ax, 'FontSize', fontsize);
legend('$\overline{STI}$ (103 SUs)', '$\overline{STI}$ (100 MUAs)', ...
    '$\overline{STI}$ (150 MUAs)', '$\overline{STI}$ (Ref)');

% title('Comparing Speech Transmission Index');




