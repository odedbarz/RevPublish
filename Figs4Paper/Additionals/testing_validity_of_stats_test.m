% testing_validity_of_stats_test.m
%
% Campare standard errors of mean for independent vs. overlapping normally
% distributed samples

% Test #1:
% 100 independent samples, each of size 30
xi = randn(30,100);
stderr_i = std(mean(xi))

% Overlapping samples: 100 subsets of size 30 from a sample of size 45
xs = randn(45,1);  % initail samples of 45
idx = zeros(30,100);
for k = 1:100,
    idx(:,k) = randperm(45,30);
end
xo = xs(idx); % 100 overlapping subsets of size 30
stderr_o = std(mean(xo))



%% Test #2:
% 100 independent samples, each of size 30
n_smp = 30;
resmp = 2;
xi = randn(n_smp,resmp);
stderr_i = std(mean(xi));

% Overlapping samples: 100 subsets of size 30 from a sample of size 45
xs = randn(n_smp*resmp,1);  % initail samples of 45
idx = zeros(n_smp,resmp);
for k = 1:resmp
    idx(:,k) = randperm(n_smp*resmp,n_smp);  % <<- stample from 30*100 instead of 45
end
xo = xs(idx); % 100 overlapping subsets of size 30
stderr_o = std(mean(xo));

fprintf('\n === Test #2: ===\n');
fprintf('Number of iid parameters: %d*%d = %d\n', n_smp, resmp, n_smp*resmp);
fprintf(' - stderr_i: %g\n', stderr_i);
fprintf(' - stderr_o: %g\n', stderr_o);

% Returns a test decision for the null hypothesis that the data in x comes 
% from a normal distribution with mean equal to zero and unknown variance, 
% using the one-sample t-test.
fprintf('\n Paired t-test:\n');
[h, p] = ttest(xs(:), xo(:));
fprintf(' - h: %d, p: %.3f\n', h, p);



%% Test #3:
% A simplified code that represent the reconstruction process with its
% resampling process (11 out of 12).
n_smp = 12;
resmp = n_smp;
zall = randn(1, n_smp);         % ZALL: one iid sample of size N_SMP
zs = zeros(1, resmp);
for k = 1:resmp
    idx = randperm(n_smp,n_smp-1);
    zs(k) = mean(zall(idx));    % resample ZALL over and over again
end

fprintf('\n === Test #3: ===\n');
fprintf(' - mean(zall): %g\n', mean(zall));
fprintf(' - mean(zs)  : %g\n', mean(zs));
fprintf(' - std(zall)   : %g\n', std(zall));
fprintf(' - std(zs)   : %g\n', std(zs));

% Returns a test decision for the null hypothesis that the data in x comes 
% from a normal distribution with mean equal to zero and unknown variance, 
% using the one-sample t-test.
fprintf('\n Paired t-test:\n');
[h, p] = ttest(zall, zs);
fprintf(' - h: %d, p: %.3f\n', h, p);




%% Test #4 -- load data
setup_environment('../../');
path_2_file = load.path_to_data('Stats');
fn = 'CCstats_TEST(1-12)_MUA_(09-Jun-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';
dummy = load(fullfile(path_2_file, fn));
stats.CCall = dummy.CCstats;
stats.units_all = dummy.CCunits;

fn = 'CCstats_TEST(1-6)_MUA_(09-Jun-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(6)_lags(30)ms_cau(0)_trainDRR(3).mat';
dummy = load(fullfile(path_2_file, fn));
stats.CC1 = dummy.CCstats;
stats.units_1 = dummy.CCunits;

fn = 'CCstats_TEST(7-12)_MUA_(09-Jun-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(6)_lags(30)ms_cau(0)_trainDRR(3).mat';
dummy = load(fullfile(path_2_file, fn));
stats.CC2 = dummy.CCstats;
stats.units_2 = dummy.CCunits;


fn = 'CCstats_TEST(1-11)_MUA_(10-Jun-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(11)_lags(30)ms_cau(0)_trainDRR(3).mat';
dummy = load(fullfile(path_2_file, fn));
stats.CC3 = dummy.CCstats;
stats.units_3 = dummy.CCunits;


% Just checking
isequal(stats.units_all, stats.units_1)
isequal(stats.units_2, stats.units_1)



%% Test #4
dry_idx = 1;
CC1 = squeeze(stats.CC1(:,dry_idx,:));
CC2 = squeeze(stats.CC2(:,dry_idx,:));
CC3 = squeeze(stats.CC3(:,dry_idx,:));
CCall = squeeze(stats.CCall(:,dry_idx,:));

figure(100);
clf;
h = aux.violinplot([std(CCall,[],2), std(CC1,[],2), std(CC2,[],2), std(CC3,[],2)],...
    {'All', 'Group 1 (1-6 utter.)', 'Group 2 (7-12 utter.)', 'Group 3 (1-11 utter.)'});
hold off
title('Comparing SD of CCs for Various Utterance Groups');
ylabel('Standard Deviation (SD)');
set(gca, 'FontSize', 32);

T = array2table([std(CCall,[],2), std(CC1,[],2), std(CC2,[],2), std(CC3,[],2)],...
    'VariableNames', {'stdCCall', 'stdCC1', 'stdCC2', 'stdCC3'});
disp(T)




