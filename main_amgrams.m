% main_amgrams.m

clc
verbose = 1;

% Load data into workspace
% analyze_setup;
data = load('D:\Dropbox (Partners HealthCare)\codeOnCloud\RevPublish\_data\Analysis\analyzed_MUA_CC_(04-Nov-2021)_unit_list(103).mat');


%% Prepare
figure(1);
t = data.spec_st.t;
f = data.spec_st.f;

Sstim_dry = data.spec_st.Sft{3};
Sstim_drr = data.spec_st.Sft{4};

dt = diff(data.spec_st.t(1:2));     % temporal sampling rate of the spectrogram

spec.plot_spectrogram(t, f, Sstim_dry);
% title(sprintf('Label: %s\\%%', data.spec_st.labels{k}(1:end-1)));



%% Create the AM-gram for one selected frequency
am_nf = 200;
fband = 26;
[am_dry, ~, am_f] = am_one_freq_band(Sstim_dry(fband,:), dt, am_nf);

%
drr_2_use = 5
Sdrr = data.Sest(:,:,drr_2_use);
am_est = am_one_freq_band(Sdrr(fband,:), dt, am_nf);

%
figure(10);
plot(am_f, [am_dry, am_est]);
% xlim([0, 0.5*1/dt]);
xlabel('Frequency (Hz)');
legend('Dry-stimulus AM', 'Est AM');



%% Calculate AMs for all frequency bands
am_nf = 200;
nf = size(Sstim_dry,1);
n_splits = data.splits.n_splits;

Sam_dry = zeros(n_splits, nf, am_nf);
Sam_drr = zeros(n_splits, nf, am_nf);
Xam_dry = zeros(n_splits, nf, am_nf);
Xam_est = zeros(n_splits, nf, am_nf);

drr_2_use = 5
Sdry = data.Sest(:,:,1);
Sdrr = data.Sest(:,:,drr_2_use);

for sp = 1:n_splits
    split_idx = sp == data.splits.idx;
    for k = 1:size(Sstim_dry,1)
        [Sam_dry(sp,k,:), ~, am_f] = am_one_freq_band(Sstim_dry(k,~split_idx), dt, am_nf);
        Sam_drr(sp,k,:)            = am_one_freq_band(Sstim_drr(k,~split_idx), dt, am_nf);
        Xam_dry(sp,k,:)            = am_one_freq_band(Sdry(k,~split_idx), dt, am_nf);
        Xam_drr(sp,k,:)            = am_one_freq_band(Sdrr(k,~split_idx), dt, am_nf);    
    end
end

Sam_dry_mean = squeeze(mean(Sam_dry,1));
Sam_dry_std  = squeeze(std (Sam_dry,1));
Sam_drr_mean = squeeze(mean(Sam_drr,1));
Sam_drr_std  = squeeze(std (Sam_drr,1));
Xam_dry_mean = squeeze(mean(Xam_dry,1));
Xam_dry_std  = squeeze(std (Xam_dry,1));
Xam_drr_mean = squeeze(mean(Xam_drr,1));
Xam_drr_std  = squeeze(std (Xam_drr,1));


%%
figure(20);
% ax = subplot(1,2,1);
% % imagesc(log10(Xam_dry));
% spec.plot_spectrogram(am_f, f, Xam_dry_mean, 'ax', gca);
% xlabel('AM Frequency (Hz)');
% colorbar;
% % xlim([1, 0.5*1/dt]);

ax = subplot(1,2,1);
% imagesc(log10(Xam_est));
spec.plot_spectrogram(am_f, f, max(0,Xam_dry_mean-Sam_dry_mean), 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
% xlim([1, 0.5*1/dt]);

ax(2) = subplot(1,2,2);
% imagesc(max(0, log10(Xam_est)-log10(Xam_dry)));
spec.plot_spectrogram(am_f, f, max(0,Xam_drr_mean-Sam_drr_mean), 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
% xlim([1, 0.5*1/dt]);

linkaxes(ax);


%%
k = 23;
xdry = Xam_dry_mean(k,:)' - Sam_dry_mean(k,:)';
xdry_std = Xam_dry_std(k,:)' - Sam_dry_std(k,:)';

xdrr = Xam_drr_mean(k,:)' - Sam_drr_mean(k,:)';
xdrr_std = Xam_drr_std(k,:)' - Sam_drr_std(k,:)';

figure(30);
plot(am_f, [xdry, xdrr]);
% errorplot([xdry, xdrr], );

h1 = plot(am_f, [xdry, xdrr]);
hold on
h2 = errorbar([am_f, am_f], [xdry, xdrr], [xdry_std, xdrr_std]);
h2(1).Color = h1(1).Color;
h2(2).Color = h1(2).Color;
plot(am_f, Sam_dry_mean(k,:)', 'k')
hold off



