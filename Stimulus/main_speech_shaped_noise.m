%
% main_speech_shaped_noise
%

clc

addpath('slamam\modSphericalHRTF');   
addpath('..');   


aux.FigSetup;

%% ====== LOAD the DRY Stimulus ======
% the path to where the WAV files will be saved 
path2save  = '.stim2Impal/';        

% The new file name to save
fn_dry = 'TIMIT_dry';

dummy = load('.\.database\TIMIT_concatenated_stimuli_Fs(100kHz)\_TIMIT_meta.mat');
Tmeta = dummy.Tmeta;
Trir  = dummy.Trir;

% Load the concatenated file
[yv, Fs] = audioread('.\.database\TIMIT_concatenated_stimuli_Fs(100kHz)\TIMIT_concatenate.wav');


%% Get the envelope (after 
foct = 1000;     % (Hz) the center frequency of the octave-band filter
[yenv, yoct] = filter_speech_env(yv, 'foct', foct);
assert(size(yenv,2) == 1, '--> ERROR: <yenv> must be a MONO signal!');


%% Plot
len_y = length(yv);
tt = linspace(0, len_y/Fs, len_y)';
ff = linspace(0, Fs, len_y)';

figure(51);
clf;
subplot(2,1,1);
plot(tt, [yoct, yenv]);
xlim([1, 1.5]);
xlabel('Time (ms)');
ylabel('Amplitude');
legend('$y_{BPF}$', '$ENV(y_{BPF})$');
subplot(2,1,2);
plot(ff, 20*log10(abs(fft([yoct, yenv]))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)'); % $20log_{10}(DFT(|\cdot|))$'
xlim([0,8e3]);


%% Speech spectrum with random phase + AM with the envelope
%{
% Keep the speech spectrum but randomize the phase
% rphs: random phase
Yv = fft(yv);

isodd = size(yv,1)/2 ~= fix(size(yv,1)/2);

if isodd    % Yv is ODD length
    dummy = Y(2:end);          % +1: don't take the DC
    Yv_half = dummy(1:end/2);
    Yv_dc     = Yv(1);
    
    % Randomize the phase
    rnd_phs    = 2*pi*rand(size(Yv_half,1),1);     % the random phase
    Yrphs_half = abs(Yv_half) .* exp(1j.*rnd_phs);

    % Stitch it all together
    Yrphs = [Yv_dc; Yrphs_half; conj(Yrphs_half(end:-1:1))];
    yrphs = ifft( Yrphs ); %, 'symmetric' );
    assert(0 == sum(imag(yrphs)), '--> ERROR: you didn''t do a proper stitching!');

else    % Yv is EVEN length
    Yv_half   = Yv(2:end/2);          % +1: don't take the DC
    Yv_dc     = Yv(1);
    Yv_middle = Yv(end/2+1);
    assert(isreal(Yv_middle), '--> ERROR: for an EVEN length signal, the center DFT coefficient is always real!');
    
    % Randomize the phase
    rnd_phs    = 2*pi*rand(size(Yv_half,1),1);     % the random phase
    Yrphs_half = abs(Yv_half) .* exp(1j.*rnd_phs);

    % Stitch it all together
    Yrphs = [Yv_dc; Yrphs_half; Yv_middle; conj(Yrphs_half(end:-1:1))];
    yrphs = ifft( Yrphs ); %, 'symmetric' );
    assert(0 == sum(imag(yrphs)), '--> ERROR: you didn''t do a proper stitching!');
    
end
%}

yrphs = shaped_speech_noise(yv);

% AM with speech envelope and fine structure of the Speech-with-random-phase
yrphs_env_mono = yrphs .* yenv;


%% Plot
figure(55);
clf;
subplot(2,1,1);
plot(tt, [yv, yrphs]);
xlim([1, 1.5]);
xlabel('Time (ms)');
ylabel('Amplitude');
legend('$y$', '$y$ random phase');
subplot(2,1,2);
plot(ff, 20*log10(abs(fft([yv, yrphs]))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)'); % $20log_{10}(DFT(|\cdot|))$'
xlim([0,8e3]);

%
figure(60);
clf;
subplot(2,1,1);
plot(tt, [yv/rms(yv), yrphs_env_mono/rms(yrphs_env_mono)]);
xlim([1, 1.5]);
xlabel('Time (ms)');
ylabel('Amplitude');
legend('$y$', '$s(t)=y_{rnd-phs}\times y_{env}$');
subplot(2,1,2);
plot(ff, 20*log10(abs(fft([yv/rms(yv), yrphs_env_mono/rms(yrphs_env_mono)]))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)'); % $20log_{10}(DFT(|\cdot|))$'
xlim([0,8e3]);



%% Convolve with the reverberation RIR
Y_rphs_env = ConvRIR(yrphs_env_mono, Trir.rir);





%% Split & SAVE the WAV files for Impale
fn_generic = sprintf('ns_Spch_fc(%gHz)', foct);
max_seg_time = 1;   % (sec)

disp(Trir);

opts.Fs         = Fs;
opts.fn_generic = fn_generic;
opts.path2save  = path2save;
opts.max_time_sec= 1.0;   % sec
opts.rms        = 'win';  

opts.params.name = {'Dist', 'Revb'};
opts.params.val = itr.Trir{:,[2,1]};
opts.params.val(:,[1,2]) = ceil(10 * opts.params.val(:,[1,2]));    % avoid the decibel numbers (1.5m --> 15)

% rms_values = split_wav_to_Impale( Y_rphs_env, Fs, Trir, inner, outer, fn_generic, path2save, max_seg_time );
rms_values = split_wav_to_Impale( Y_rphs_env, Trir, opts );

fprintf('--> rms_values:\n');
disp(rms_values);



%% Save all meta information about the stimuli (the "meta" data)
meta_fn = [path2save, '_', fn_generic, '_meta', '.mat'];
save( meta_fn, 'Fs', 'foct', 'fn_dry' );
aux.cprintf('cyan', '--> META-data file was saved (%s)\n', meta_fn);







