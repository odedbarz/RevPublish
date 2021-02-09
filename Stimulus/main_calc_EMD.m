%
% main_calc_EMD.m
%

clc

% Get all WAV files in ONE table (tbl_wav)
main_find_TIMITs_in_dry;   % --> tbl_wav

% Calculate EMD
% %{
fprintf('\n');
fprintf('--> Starting EMD ...\n');

% Add the IMF column
tbl_imf = [tbl_wav, array2table( cell(n_rows,1), 'VariableNames', {'imf'} )];

% Nstd    = 0.2*rms(y);
NR      = 10;       %500;
MaxIter = 1000;     %5000;
SNRFlag = 1;

for k = 1:size(tbl_imf,1)
    yk = tbl_imf.wav_short{k};
    Nstd    = 0.2*rms(yk);

    [imf, ~] = ceemdan.ceemdan_v2014(yk, Nstd, NR, MaxIter, SNRFlag);
    tbl_imf.imf{k} = imf';
    
    
end
%}

% % Load the IMFs instead of (calculating them)
% load('D:\GoogleDrive\codeOnCloud\Reverberation\Stimulus\.database\Project1\Spch_(36)sec\tbl_wav.mat');

'NOT SAVING!'
% save('tbl_imf.mat', 'tbl_imf', 'NR', 'MaxIter', 'SNRFlag');

