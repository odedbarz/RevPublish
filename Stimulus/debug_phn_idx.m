%
% debug_phn_idx.m
% 
% A quick and nasty "correction" script to adjust the indices ofthe
% phonetic enries.
% 

clc

%%
dummy = load('.\.database\TIMIT_dry\TIMIT_dry_meta.mat');
Tmeta = dummy.Tmeta;


%% Load the concatenated file
[ya, fy] = audioread(['.\.database\TIMIT_dry\TIMIT_dry.wav']);


%% Load the TIMIT file 
% For the current selected speech sample
% %{
wav_file_list{1} = 'TIMIT_train_dr1_MRWS0_SX292';
wav_file_list{2} = 'TIMIT_train_dr6_MSDS0_SI1707';
wav_file_list{3} = 'TIMIT_train_dr6_MSJK0_SX156';
wav_file_list{4} = 'TIMIT_train_dr6_FAPB0_SX163';
wav_file_list{5} = 'TIMIT_train_dr3_MWGR0_SA1';
wav_file_list{6} = 'TIMIT_train_dr1_FCJF0_SA1';
wav_file_list{7} = 'TIMIT_train_dr2_FAEM0_SA2';
wav_file_list{8} = 'TIMIT_test_dr1_FAKS0_SX313';
wav_file_list{9} = 'TIMIT_test_dr1_MDAB0_SX409';
wav_file_list{10} = 'TIMIT_test_dr7_FCAU0_SX407';
wav_file_list{11} = 'TIMIT_test_dr6_FDRW0_SX293';
wav_file_list{12} = 'TIMIT_train_dr7_MFXV0_SX375';    
%}

path_abs = 'D:\DataBases\TIMIT\';
fn_list = {'SX292', 'SI1707', 'SX156', 'SX163', 'SA1', 'SA1', 'SA2',...
    'SX313', 'SX409', 'SX407', 'SX293', 'SX375'};

for kk = 1:length(fn_list)
    
    fn = regexprep(wav_file_list{kk}, '_', '\');
    fn = [path_abs, fn];
    
    % load this file
    [x, fx] = audioread([fn, '.wav']);

    % Read the PHN file
    [phn_start, phn_end, phn] = textread([fn, '.phn'], '%d %d %s'); 

    % Read the WRD file
    [wrd_start, wrd_end, wrd] = textread([fn, '.wrd'], '%d %d %s'); 


    %% Resample 
    % Resample the signal
    [num, den] = rat(fy/fx);
    xi = resample(x, num, den);
    xi = xi/std(xi);
    fxi = fy;
    len_xi = length(xi);

    % Resample the indices
    % 1+: indices in TIMIT start at 0; convert to MATLAB notation
    phn_start_i = fix(1+ phn_start*num/den);
    phn_end_i   = fix(1+ phn_end*num/den);
    wrd_start_i = fix(1+ wrd_start*num/den);
    wrd_end_i   = fix(1+ wrd_end*num/den);


    %% 
    % yk = ya(1:len_xi);
    yk = ya; %(phn_start_i(1):phn_end_i(end));
    yk = yk/std(yk);
    [C, lags] = xcorr(yk, xi);

    [~,idx] = max(C);
    lag_xi = lags(idx);

    %{
    plot(yk, '.--', 'LineWidth', 2);
    hold on
    plot(lag_xi + (1:length(xi)), xi, 'LineWidth', 1);
    hold off
    grid on
    title(sprintf('kk: %d', kk));
    drawnow;
    %}
    
    %% Remove the silence the h# and 
    phn_start_i = phn_start_i(2:(end-1));
    phn_end_i   = phn_end_i(2:(end-1));

    %% Remove the lag
    phn_start_i = phn_start_i + lag_xi;
    phn_end_i   = phn_end_i   + lag_xi;

    % Update the phonetic indices
    Tmeta.phn_start{kk} = phn_start_i;
    Tmeta.phn_end{kk}   = phn_end_i;
    Tmeta.phn{kk}       = Tmeta.phn{kk}(2:(end-1));
    
    % Update the words indices
    Tmeta.wrd_start{kk} = wrd_start_i + lag_xi;
    Tmeta.wrd_end{kk}   = wrd_end_i   + lag_xi;
    
end


Tmeta



