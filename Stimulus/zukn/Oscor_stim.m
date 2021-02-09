function files = Oscor_stim(mod_rates, num_tokens, varargin)
% Creates a series of oscillating correlation stimuli for use in Impale

% Starting variables
t_end = 5;
Fs = 100000;
start_tok = randi(100);
phase = pi/2; % the phase between the correlated and uncorrelated noises when summed
    % in one ear
path = 'V:\ITD\Oscor\';
save_flag = 0;
play_flag = 1;
am_depth = 0;
am_freq = 4;

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

t = 0:1/Fs:t_end-1/Fs;
files = {};
for m = 1:length(mod_rates),
    for n = 1:num_tokens,
        rng(start_tok+n-1);
        % Create two uncorrelated random Gaussian noises
        A = randn(length(t),1);
        B = randn(length(t),1);
        % The amplitudes of the noises, when summed in 
        A_mot = -cos(2*pi*mod_rates(m)*t)'; % -cos so that it starts anti-correlated
            % akin to the choice of the DynITD stimuli
        B_mot = -cos(2*pi*mod_rates(m)*t-phase)'; % lags A_mot by phase
        % Sum these waveforms for one ear
        if ~mod_rates(m),
            IP_chan = A; % if the modulation rate is 0, use correlated noise
        else
            IP_chan = A.*A_mot + B.*B_mot;
        end
        CO_chan = A;
        % Apply amplitude modulation
        am = 10^(am_depth/20);
        env = 1/2*(1-am*cos(2*pi*am_freq*t))';
        IP_chan = IP_chan.*env;
        CO_chan = CO_chan.*env;
        % Normalize channels by maximum value of both
        mu = max(max(abs([IP_chan CO_chan])));
        IP_chan = IP_chan/mu;
        CO_chan = CO_chan/mu;
        % Save files
        if save_flag,
            IP_file = [path 'Oscor_IP_' ...
                '-v' num2str(mod_rates(m)) ...
                '-Rep' num2str(n)];
            CO_file = [path 'Oscor_CO_' ...
                '-v' num2str(mod_rates(m)) ...
                '-Rep' num2str(n)];
                EPLwavwrite(IP_chan,Fs,16,IP_file,'rms_dB',db(sqrt(mean(IP_chan.^2))));
            disp(IP_file);
            EPLwavwrite(CO_chan,Fs,16,CO_file,'rms_dB',db(sqrt(mean(CO_chan.^2))));
            disp(CO_file);
            files = [files {IP_file, CO_file}];
        end
        if play_flag,
            sound([IP_chan CO_chan],Fs);
            files = [IP_chan CO_chan];
        end
    end
end
    