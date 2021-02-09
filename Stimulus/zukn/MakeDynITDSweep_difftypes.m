function [files] = MakeDynITDSweep(mod_rates, num_tokens, amp_ITD, varargin)
% Run MotSource_ITDSweep multiple times to create a set of stimuli (saved
% with EPLwavread) for the set of mod_rates and num_tokens specified, with
% the ITD sweep range of +/-amp_ITD

% Initial variables
start_tok = randi(100);
Fs = 100000;
path = 'V:\ITD\DynITD\sweeps\';

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

%%% TODO (7/28/14): Setup different sweep types, saving those, also change filenames
files = {}; % Create cell array to store filenames
for m = 1:length(mod_rates),
    for tok = 1:num_tokens,
        % Compute stimulus duration
        sweep_dur = 1/(2*mod_rates(m));
        % Create stimuli
        %% Constant sweeps
        s_pos = MotSource_ITDSweep(0, amp_ITD, mod_rates(m), 'tok', start_tok+tok-1, 'sweep_type', 'constant', 'Fs', Fs);
        s_neg = MotSource_ITDSweep(0, -amp_ITD, mod_rates(m), 'tok', start_tok+tok-1, 'sweep_type', 'constant', 'Fs', Fs);
        % Compute rms
%         rms_CO_pos = sqrt(mean(s_pos(:,2).^2));
%         rms_IP_pos = sqrt(mean(s_pos(:,1).^2));
%         rms_CO_neg = sqrt(mean(s_neg(:,2).^2));
%         rms_IP_neg = sqrt(mean(s_neg(:,1).^2));
        % Compute rms, not including silence at end of stim
        rms_CO_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),1).^2));
        rms_CO_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),1).^2));
        % Name stimulus
        IP_chan_pos = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(3) ...
            '_-Rep' num2str(tok)];
        CO_chan_pos = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(3) ...
            '_-Rep' num2str(tok)];
        IP_chan_neg = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-3) ...
            '_-Rep' num2str(tok)];
        CO_chan_neg = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-3) ...
            '_-Rep' num2str(tok)];
        % Save the stimulus file
        EPLwavwrite(s_pos(:,1),Fs,16,IP_chan_pos,'rms_dB',db(rms_IP_pos));
        EPLwavwrite(s_pos(:,2),Fs,16,CO_chan_pos,'rms_dB',db(rms_CO_pos));
        EPLwavwrite(s_neg(:,1),Fs,16,IP_chan_neg,'rms_dB',db(rms_IP_neg));
        EPLwavwrite(s_neg(:,2),Fs,16,CO_chan_neg,'rms_dB',db(rms_CO_neg));

        files = [files, {IP_chan_pos}, {CO_chan_pos}, {IP_chan_neg}, {CO_chan_neg}];
        disp([IP_chan_pos '; ' CO_chan_pos]);
        disp([IP_chan_neg '; ' CO_chan_neg]);
        
        %% Steps
        s_pos = MotSource_ITDSweep(0, amp_ITD, mod_rates(m), 'tok', start_tok+tok-1, 'sweep_type', 'step','Fs',Fs);
        s_neg = MotSource_ITDSweep(0, -amp_ITD, mod_rates(m), 'tok', start_tok+tok-1, 'sweep_type', 'step','Fs',Fs);
%         % Compute rms
%         rms_CO_pos = sqrt(mean(s_pos(:,2).^2));
%         rms_IP_pos = sqrt(mean(s_pos(:,1).^2));
%         rms_CO_neg = sqrt(mean(s_neg(:,2).^2));
%         rms_IP_neg = sqrt(mean(s_neg(:,1).^2));
        % Compute rms, not including silence at end of stim
        rms_CO_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),1).^2));
        rms_CO_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),1).^2));
        % Name stimulus
        IP_chan_pos = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(2) ...
            '_-Rep' num2str(tok)];
        CO_chan_pos = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(2) ...
            '_-Rep' num2str(tok)];
        IP_chan_neg = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-2) ...
            '_-Rep' num2str(tok)];
        CO_chan_neg = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-2) ...
            '_-Rep' num2str(tok)];
        % Save the stimulus file
        EPLwavwrite(s_pos(:,1),Fs,16,IP_chan_pos,'rms_dB',db(rms_IP_pos));
        EPLwavwrite(s_pos(:,2),Fs,16,CO_chan_pos,'rms_dB',db(rms_CO_pos));
        EPLwavwrite(s_neg(:,1),Fs,16,IP_chan_neg,'rms_dB',db(rms_IP_neg));
        EPLwavwrite(s_neg(:,2),Fs,16,CO_chan_neg,'rms_dB',db(rms_CO_neg));

        files = [files, {IP_chan_pos}, {CO_chan_pos}, {IP_chan_neg}, {CO_chan_neg}];
        disp([IP_chan_pos '; ' CO_chan_pos]);
        disp([IP_chan_neg '; ' CO_chan_neg]);
        
        %% Static
        dur_stat = 0.5*1/(2*mod_rates(m));
        % Create static stimuli by specifying the edges of the ITD range as
        % the "center ITD", and make sweep without modulation for half the
        % duration of a typical sweep (dur_stat)
        s_pos = MotSource_ITDSweep(-amp_ITD, 0, 0, 'tok', start_tok+tok-1, 't_end', dur_stat,'Fs',Fs);
        s_neg = MotSource_ITDSweep(amp_ITD, 0, 0, 'tok', start_tok+tok-1, 't_end', dur_stat,'Fs',Fs);
%          % Compute rms
%         rms_CO_pos = sqrt(mean(s_pos(:,2).^2));
%         rms_IP_pos = sqrt(mean(s_pos(:,1).^2));
%         rms_CO_neg = sqrt(mean(s_neg(:,2).^2));
%         rms_IP_neg = sqrt(mean(s_neg(:,1).^2));
        % Compute rms, not including silence at end of stim
        rms_CO_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_pos = sqrt(mean(s_pos(1:floor(sweep_dur*Fs),1).^2));
        rms_CO_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),2).^2));
        rms_IP_neg = sqrt(mean(s_neg(1:floor(sweep_dur*Fs),1).^2));
        % Name stimulus
        IP_chan_pos = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(1) ...
            '_-Rep' num2str(tok)];
        CO_chan_pos = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(1) ...
            '_-Rep' num2str(tok)];
        IP_chan_neg = [path 'DynITDSweep_IP' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-1) ...
            '_-Rep' num2str(tok)];
        CO_chan_neg = [path 'DynITDSweep_CO' ...
            '_-v' num2str(mod_rates(m)) ...
            '_-sw' num2str(-1) ...
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