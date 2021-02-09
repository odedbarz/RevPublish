function [files] = MakeDynITDStimuli_interp(mod_rates, cntr_ITD, num_tokens, amp_ITD, varargin)
% Run MotSource_interp multiple times and save each stimulus waveform using
% EPLwavread.  These create dynamic ITD stimuli without AM.
% amp_ITD should be positive

% Initial variables
start_tok = 42;
Fs = 100000;
t_end = 5;
path = 'V:\ITD\DynITD\noise_interp_alwaysdecrdelay\'; %%% EDIT (11/5/14)

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

files = {}; % Create cell array to store filenames
for m = 1:length(mod_rates),
    for c = 1:length(cntr_ITD),
        for tok = 1:num_tokens,
            % Create stimulus
            %%% (11/5/14) Filter before and after ITD modulation.
            %%% (11/5/14) Vary the delays in both channels so the delay is
            %%% always decreasing
            s = MotSource_interp(cntr_ITD(c),amp_ITD,mod_rates(m),2,-Inf,0,0,...
                't_end',t_end,'motion_type','triangle','Fs',Fs,'tok',start_tok+tok-1,'both_delay',1,'flt_noise',1);
            % Compute rms
            rms_CO = sqrt(mean(s(:,2).^2));
            rms_IP = sqrt(mean(s(:,1).^2));
            % Name stimulus
            IP_chan = [path 'DynITD_IP' ...
                '_-v' num2str(mod_rates(m)) ...
                '_-c' num2str(cntr_ITD(c)) ...
                '_-Rep' num2str(tok)];
            CO_chan = [path 'DynITD_CO' ...
                '_-v' num2str(mod_rates(m)) ...
                '_-c' num2str(cntr_ITD(c)) ...
                '_-Rep' num2str(tok)];
            % Save the stimulus file
            EPLwavwrite(s(:,1),Fs,16,IP_chan,'rms_dB',db(rms_IP));
            EPLwavwrite(s(:,2),Fs,16,CO_chan,'rms_dB',db(rms_CO));
            
            files = [files, {IP_chan}, {CO_chan}];
            disp([IP_chan '; ' CO_chan]);
        end
    end
end