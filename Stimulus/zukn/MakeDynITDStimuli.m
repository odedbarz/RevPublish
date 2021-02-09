function [files] = MakeDynITDStimuli(mod_rates, cntr_ITD, num_tokens, amp_ITD, varargin)
% Run MotSource_ITDnAM multiple times and save each stimulus waveform using
% EPLwavread.  These create dynamic ITD stimuli without AM.
% amp_ITD should be positive

% Initial variables
start_tok = randi(100);
Fs = 100000;
t_end = 5;
path = 'V:\ITD\DynITD\mod_examine_alwaysdecrdelay\';

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
            s = MotSource_ITDnAM(cntr_ITD(c),amp_ITD,mod_rates(m),2,-Inf,0,0,...
                't_end',t_end,'motion_type','triangle','Fs',Fs,'tok',start_tok+tok-1);
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