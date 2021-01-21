% load_py.m
%
%

setup_environment('../');

mat_fn = 'analyzed_stimulus.mat';
mat_full_fn = fullfile( load.path_to_data('data'), 'Analysis', mat_fn);
data = load( mat_full_fn );
sr = (data.sr);
duration_seconds = double(data.duration_seconds);
fmin = double(data.fmin);
fmax = double(data.fmax);
n_fft = (data.n_fft);
t = data.t;
y = double(data.y);
t_ = data.t_;
F0 = data.F0;
harm = data.harm;
voiced_flag = data.voiced_flag;
hop_length = double(data.hop_length);



%%
figure(1);
clf;
plot(t, 500*y, t_, F0)



nf_ = size(harm,1);    % # of frequencies
nt_ = size(harm,2);    % # of time steps

f_ = sr/2 * (0:nf_-1)'/nf_;
% f_ = 8192 * (0:nf_-1)'/nf_;
t_ = duration_seconds * (0:nt_-1)'/nt_;

% harm = nan(size(harm));
harm(:, 1 == voiced_flag) = harm(:, 1 == voiced_flag);

    
figure(2);
clf;
hz = 1;       % set the units: 1 (Hz) or 1e-3 (kHz)
imgh = imagesc(t_, hz*f_, harm);
hold on
% plot(t_, hz*F0-fmin, 'w')
plot(t_, hz*F0, 'w')
hold off
set(gca,'YDir','normal')
ylim(hz*[fmin, fmax]);



%%
window = n_fft;
[S, f, t] = spectrogram(y, window, hop_length, n_fft, sr);
Sdb = db(abs(S));
Sdb = max(0, Sdb + 0);

figure(3);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)
imagesc(t, hz*f, Sdb);
colorbar;
set(gca,'YDir','normal')
ylim(hz*[0, 1e3]);
% set(gca, 'YScale', 'log')

hold on
% plot(t_, hz*F0-fmin, 'w')
plot(t_, hz*F0, 'w')
hold off









