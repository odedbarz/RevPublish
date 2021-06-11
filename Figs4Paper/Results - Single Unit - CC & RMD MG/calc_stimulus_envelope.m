function [Yenv, env_st] = calc_stimulus_envelope(Y, CF, fs_psth, fs_stim, win_size_ms)
%
%
% function [Yenv, oct] = calc_stimulus_envelope(Y, CF, fs_psth, fs_stim, win_size_ms)
%
%

if 5 > nargin
    win_size_ms = 25;    
end

% Envelope upsampling ratio
env_upsmp = 2; %4;

% Octave BPF order
Noct = 3;

%% ENVELOPEs of the STIMULUS (as a function of the neuron's CF)
if isnan(CF) || (0 == CF)   %|| (CF>8e3)  % is there a valid CF?
    % ##############################
    % ##### Broadband Envelope #####
    % ##############################
    env_st.win_size_ms = win_size_ms;   % (ms)
    env_st.n_win = ceil( 1e-3*env_st.win_size_ms * fs_stim);  % # of samples for the window
    env_st.win_type = 'hanning';
    lpf_win = eval(sprintf('%s(%d)', env_st.win_type, env_st.n_win));

    % Rectify the matrix with all the stimuli
    Y_rectify = max(0, Y);

    % Create the envelope; filter the stimuli with the window
    Yenv = filtfilt(lpf_win, 1, Y_rectify);

    % Resampling ratio
    [num, den] = rat(fs_psth/fs_stim);


else    % 0 < CF < 8kHz
    % ################################
    % ##### Octave-BPF + Hilbert #####
    % ################################

    % Upsample the stimulus for the octave BPF
    Yup = resample(Y, env_upsmp, 1);
    stim_fs_up = env_upsmp * fs_stim;

    % Create the octave BPF
    try
        [env_st.B, env_st.A] = octave.octdsgn(CF, stim_fs_up, Noct);
    catch EM
        fprintf('\n--> Error Message: %s\n', EM.message);
        fprintf('\n--> There was an error in the octave BPF at\n');
        fprintf('---> CF     : %g Hz\n', 1e-3*CF);
        error('');
    end
    Ybpf = filter(env_st.B, env_st.A, Yup);

    % Create the envelope
    Yenv = envelope( Ybpf );

    % Resampling ratio
    [num, den] = rat(fs_psth/stim_fs_up);

    env_st.resampling.num = num;
    env_st.resampling.den = den;
end


% Resample (downsample) the envelope to the PSTH's sampling rate
Yenv = resample(Yenv, num, den);
env_st.Yenv = Yenv;


if 1 < nargout
    Yenv = env_st.Yenv;
end






