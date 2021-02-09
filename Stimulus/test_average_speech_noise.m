%
% test_average_speech_noise.m
%
% Test the function that generates the averaged-spectrum random-phase stimulus.
%

clc

files.path   = {'D:\DataBases\TIMIT\Female', 'D:\DataBases\TIMIT\Males'}; 
files.all = extract_wav_files( files.path );

[Ywav, fs] = load_wav_files( files.all );
fprintf('--> fs: %g Hz\n', fs);

[yavg, num_wav, Yfft] = averaged_speech_noise( Ywav );


%% Plot the results
len_y = size(Yfft,1);
tt = linspace(0, 1/fs*len_y, len_y)';
ff = linspace(0, fs, len_y)';

figure(11);
clf;
subplot(2,1,1);
plot(ff, 20*log10( abs(Yfft) ));
title(sprintf('A batch of %d spectrum responses', num_wav),...
    'Interpreter', 'latex', 'FontSize', 24);
ylabel('Spectrum (dB)', 'Interpreter', 'latex');
xlabel('Samples', 'Interpreter', 'latex');
hold on
plth = plot(ff, 20*log10( abs(fft(yavg)) ), '--k');
hold off
legend(plth, '20*log10(|DFT(y)|): Avg. Spectrum Response');

subplot(2,1,2);
plot(tt, yavg );
xlabel('Time (sec)', 'Interpreter', 'latex');
ylabel('Amp.', 'Interpreter', 'latex');
title('$y_{avg}$', 'Interpreter', 'latex', 'FontSize', 24);






