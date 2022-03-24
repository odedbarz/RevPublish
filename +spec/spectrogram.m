function [Sft, spec_st, f, t] = spectrogram(y, fs, varargin)
%
%   function [sgram_st, Sft, f, t]  = spectrogram(y, fs_time, [struct] or [cell array], )
%

%% Set the inputs
p = inputParser;

addRequired(p, 'y', @isnumeric);            
addRequired(p, 'fs', @isnumeric);      % (Hz, 1x1) sampling rate of the input stimulus <y>

addOptional(p, 'n_bands', 60, @isnumeric);       
addOptional(p, 'lowfreq', 250, @isnumeric);             
addOptional(p, 'highfreq', 8000, @isnumeric);      
addOptional(p, 'overlap_ratio', 0.80, @isnumeric);      
addOptional(p, 'nw', 2, @isnumeric);                % (1x1) for the multitaper spectrogram 
addOptional(p, 'win_size_ms', 25, @isnumeric);     
addOptional(p, 'binwidth', 1, @isnumeric);     
addOptional(p, 'duration_ms', [], @isnumeric);             
addOptional(p, 'db_floor', -80, @isnumeric);             
addOptional(p, 'f_scale', 'log', @isstr);           % {'lin', 'log', 'erb'}
addOptional(p, 'amp_to_db', true, @(x) isnumeric(x) | islogical(x));      % 
addOptional(p, 'method', 'matlab', @isstr);         % {'matlab', 'tfspec', 'multitaper', *'meddis'}          
addOptional(p, 'rabbit_flag', false, @isstr);           
addOptional(p, 'apply_sync_filter', false, @(x) islogical(x) | isnumeric(x));    
addOptional(p, 'calc_envelope', true, @(x) islogical(x) | isnumeric(x));    
addOptional(p, 'fignum', [], @isnumeric);           % (1x1) if not empty, plot stuff (for debug)

parse(p, y, fs, varargin{:});

spec_st    = p.Results;
spec_st.fs = fs;            % (Hz)
fignum     = p.Results.fignum;      



%% Frequency scale (linear of logarithmic scaled)
if strcmpi('log', spec_st.f_scale)
    spec_st.f = logspace(log10(spec_st.lowfreq), log10(spec_st.highfreq), spec_st.n_bands)';
    spec_st.f(1) = round(spec_st.f(1));
    spec_st.f(end) = round(spec_st.f(end));

elseif strcmpi('lin', spec_st.f_scale)
    spec_st.f = linspace(spec_st.lowfreq, spec_st.highfreq, spec_st.n_bands)'; 

elseif strcmpi('erb', spec_st.f_scale)    
    % The CF filters in the cochlea
    spec_st.f = ERBSpace(spec_st.lowfreq, spec_st.highfreq, spec_st.n_bands, spec_st.rabbit_flag);     % (also implemented in MakeERBFilters.m)
    spec_st.f = spec_st.f(end:-1:1);
    spec_st.f(1) = round(spec_st.f(1));
    spec_st.f(end) = round(spec_st.f(end));
    
else
    error('--> ERROR: unrecognized F_SCALE (''%s'')!!\n', spec_st.f_scale);
end


%% Calculate the window
fs_step = 1/units.ms2sec(spec_st.binwidth);


%% SPECTROGRAM
% Select the spectrogram method
switch lower(spec_st.method)
    case 'matlab'
        % Sets the tapered window size such that the output spectrum Sx has
        % the desired time step of BINWIDTH
        %spec_st.n_win   = 5 * spec_st.fs / fs_step;
        assert( ~isempty(spec_st.win_size_ms) && ~isnan(spec_st.win_size_ms) );        
        spec_st.n_win = fix(1e-3*spec_st.win_size_ms * spec_st.fs);
        
        % Overwrite the WIN_SIZE_MS        
        %spec_st.win_size_ms = units.sec2ms(spec_st.n_win/spec_st.fs);       % (ms)
        spec_st.overlap = ceil( spec_st.overlap_ratio * spec_st.n_win ); 
        spec_st.win     = hanning( spec_st.n_win );

        % Pad with zeros by one extra window-size to get the desired duration (samples)
        y = [y(:); zeros(spec_st.n_win-1,1)];
        
        % '### [sgram.m]: using MATLAB''s spectrogram ###'
        [Sx, spec_st.f, spec_st.t] = spectrogram(y, spec_st.win, spec_st.overlap, spec_st.f, spec_st.fs);
    
    case 'stft'
        spec_st.n_win   = round( spec_st.win_size_ms * spec_st.fs / fs_step );
        spec_st.overlap = ceil( spec_st.overlap_ratio * spec_st.n_win ); 
        spec_st.win     = hanning( spec_st.n_win );
        
        %'### [sgram.m]: using MY spectrogram function! ###'       
        [Sx, spec_st.f, spec_st.t] = spec.stft(y, spec_st.fs,...
            'f', spec_st.f, ...
            'abs_spectrum', 0,...
            'win', spec_st.win, ...
            'binwidth', spec_st.binwidth, ...
            'do_log', 0,...
            'remove_mean', 0 ...
        );
        %assert(isequal(spec_st.f,f));
        
    case {'gammatone', 'meddis', 'carney'}
        assert(strcmpi('erb', spec_st.f_scale), 'For gammatone filters set scale to log scale!')
        
        fs_new = 1/(1e-3*spec_st.binwidth);     % (Hz)
        [N, downsample] = rat(fs_new/spec_st.fs);
        assert(1 == N);
        
        [Sx, spec_st.f, spec_st.t] = Stim2ANF( y,...
            'method', spec_st.method, ...
            'calc_envelope', spec_st.calc_envelope, ...
            'Fs', spec_st.fs, ...
            'downsmp', downsample, ...
            'Nch', spec_st.n_bands, ...
            'apply_sync_filter', p.Results.apply_sync_filter,...
            'lowfreq', spec_st.lowfreq, ...
            'highfreq', spec_st.highfreq  ...
        );    
                
        
    case 'multitaper'
        spec_st.n_win   = round( spec_st.win_size_ms * spec_st.fs / fs_step );
        
        '### [sgram.m]: using MULTITAPER cochleogram!! ###'
        [Sx, spec_st.f, spec_st.t] = spec.multitaper(y, spec_st.fs,...
            'f', spec_st.f, ...
            'binwidth', spec_st.binwidth, ...
            'n_win', spec_st.n_win, ...
            'nw', spec_st.nw, ...
            'do_log', 0,...
            'remove_mean', 0 ...
        );
        %assert(isequal(spec_st.f,f));
        
    otherwise
        error('--> ERROR at [sgram.m]: unrecognized METHOD string (method: %s)!!!', lower(method));
end



%% Make sure that the spectrogram has the right temporal length
binwidth_ = 1e3* diff(spec_st.t(1:2));       % (ms)
if norm(binwidth_ - spec_st.binwidth) > 100*eps(class(binwidth_))
    warning('--> [spec.spectrogram]: desired time-axis frequency does NOT equal the actual one!!');
    spec_st.binwidth = binwidth_;
end

if isempty(spec_st.duration_ms)
    % If not given, calculate the duration time. If given, complete (or
    % truncate) to the desired duration time
    spec_st.duration_ms = size(Sx, 2)*spec_st.binwidth;
end


spec_st.t = (0:(size(Sx,2)-1)) * (1e-3*spec_st.binwidth);
spec_st.n_time = length(spec_st.t);           % (samples)     






%% Normalize for the dB scale
if spec_st.amp_to_db
    % Get the absolute value of the spectrogram
    Sx_abs = abs(Sx);

    spec_st.max_Sft = max(Sx_abs(:));

    % Apply a threshold to the spectrogram
    spec_st.Sft = 20*log10( (eps + Sx_abs)/(eps + spec_st.max_Sft) ) - spec_st.db_floor;
    spec_st.thr = 0;
    spec_st.Sft = max(spec_st.thr, spec_st.Sft);
    Sft         = spec_st.Sft;
else
    Sft = max(0, Sx);
end


if 1 < nargout
    t   = spec_st.t;
    f   = spec_st.f;
end




%% Plot
if ~isempty(fignum)
    spec.plot_spectrogram(spec_st.t, spec_st.f, spec_st.Sft, fignum);
end





























