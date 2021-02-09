function [files] = MakeDynITDSweep(mod_rates, cntr_ITDs, num_tokens, amp_ITD, varargin)
% Run ITDsweepstim multiple times to create a set of stimuli (saved
% with EPLwavread) for the set of mod_rates and num_tokens specified, with
% the ITD sweep range of +/-amp_ITD

% Initial variables
start_tok = randi(100);
Fs = 100000;
path = 'V:\ITD\DynITD\sweeps_btwnuncorr\';
uncdur = 250; % duration of starting and ending uncorrelated noise
t_end = 5; % total duration of the stimulus

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

files = {}; % Create cell array to store filenames
for m = 1:length(mod_rates),
    for c = 1:length(cntr_ITDs),
        for tok = 1:num_tokens,
            % Compute stimulus duration
            sweep_dur = 1/(2*mod_rates(m));
            % Create stimuli
            %% Constant sweeps
            s_pos = ITDsweepstim(sweep_dur*1000, amp_ITD, cntr_ITDs(c), 1, uncdur, 'tok', start_tok+tok-1, ...
                'Fs', Fs, 't_end', t_end);
            s_neg = ITDsweepstim(sweep_dur*1000, amp_ITD, cntr_ITDs(c), -1, uncdur, 'tok', start_tok+tok-1, ...
                'Fs', Fs, 't_end', t_end);
            % Normalize so the maximum value of the sound is 1 or -1
            s_pos = s_pos/max(max(abs(s_pos)));
            s_neg = s_neg/max(max(abs(s_neg)));
            % Compute rms, not including silence at end of stim
            rms_CO_pos = sqrt(mean(s_pos(:,2).^2));
            rms_IP_pos = sqrt(mean(s_pos(:,1).^2));
            rms_CO_neg = sqrt(mean(s_neg(:,2).^2));
            rms_IP_neg = sqrt(mean(s_neg(:,1).^2));
            % Name stimulus
            IP_chan_pos = [path 'DynITDSweep_IP' ...
                '_-v' num2str(mod_rates(m)) ...
                '_-c' num2str(cntr_ITDs(c)) ...
                '_-Rep' num2str(tok)];
            CO_chan_pos = [path 'DynITDSweep_CO' ...
                '_-v' num2str(mod_rates(m)) ...
                '_-c' num2str(cntr_ITDs(c)) ...
                '_-Rep' num2str(tok)];
            IP_chan_neg = [path 'DynITDSweep_IP' ...
                '_-v' num2str(-mod_rates(m)) ...
                '_-c' num2str(cntr_ITDs(c)) ...
                '_-Rep' num2str(tok)];
            CO_chan_neg = [path 'DynITDSweep_CO' ...
                '_-v' num2str(-mod_rates(m)) ...
                '_-c' num2str(cntr_ITDs(c)) ...
                '_-Rep' num2str(tok)];
            % Save the stimulus file
            EPLwavwrite(s_pos(:,1),Fs,16,IP_chan_pos,'rms_dB',db(rms_IP_pos));
            EPLwavwrite(s_pos(:,2),Fs,16,CO_chan_pos,'rms_dB',db(rms_CO_pos));
            EPLwavwrite(s_neg(:,1),Fs,16,IP_chan_neg,'rms_dB',db(rms_IP_neg));
            EPLwavwrite(s_neg(:,2),Fs,16,CO_chan_neg,'rms_dB',db(rms_CO_neg));

            files = [files, {IP_chan_pos}, {CO_chan_pos}, {IP_chan_neg}, {CO_chan_neg}];
            disp([IP_chan_pos '; ' CO_chan_pos]);
            disp([IP_chan_neg '; ' CO_chan_neg]);
        end
    end
end