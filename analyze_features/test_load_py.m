% test_load_py.m
%
% Description:
% Loads analyzed data of the dry stimulus that was processed with LIBROSA
% (python).

setup_environment('../');


% Loads the MAT file (librosa)
mat_fn = 'analyzed_librosa.mat';
mat_full_fn = fullfile( load.path_to_data('data'), 'Analysis', mat_fn);
data = load( mat_full_fn );

% Extract parameters
spec = data.spec;
p = data.p;

sr = (data.stim.sr);
duration_seconds = double(data.stim.duration_seconds);
fmin = double(data.p.fmin);
fmax = double(data.p.fmax);
n_fft = (data.spec.n_fft);
t = data.stim.t;
y = double(data.stim.y);
t_ = data.p.t_;
F0 = data.p.F0;
harm = data.harm;
voiced_flag = data.p.voiced_flag;
hop_length = double(data.p.hop_length);



%%
figure(3);
clf;
hz = 1e-3;       % set the units: 1 (Hz) or 1e-3 (kHz)

img = imagesc(spec.t, hz*spec.f, spec.Sdb);
set(gca, 'YDir', 'normal');

colorbar;
ylim(hz*[0, 1e3]);
% set(gca, 'YScale', 'log')

hold on
% plot(t_, hz*F0-fmin, 'w')
plot(p.t_, hz*p.F0, 'w')
hold off













