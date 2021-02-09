function [rms_values, file_list] = split_wav_to_Impale( Y_list, opts)
%
% [rms_values, file_list] = split_wav_to_Impale( Y_list, Trir, opts)
%
% Options:
% opts.Fs         = Fs;
% opts.fn_generic = fn_generic;
% opts.path2save  = path2save;
% opts.duration_sec= 3;   % sec
% 
% opts.params.name = {'Dist', 'Revb'};
% opts.params.table_cols = [1, 2];     % [1, 2] <--> [wtypes | dist ]
%
%

n_signals = size(Y_list,1);

% % Set the input variables
% assert(size(Y_list,1) == n_signals);

% % outer & inner MUST be columns in Trir table!
% assert(1 == prod(size(Trir,2) >= outer_inner_cols),...
%     '--> ERROR in [split_wav_to_Impale.m]: Outer & inner MUST be columns in Trir table!');  
% 

assert( isfield(opts, 'Fs'), '--> ERROR in [split_wav_to_Impale.m]: OPTS must contain <Fs>!' );
assert( isfield(opts, 'fn_generic'), '--> ERROR in [split_wav_to_Impale.m]: OPTS must contain <fn_generic>!' );

if ~isfield(opts, 'path2save')
    opts.path2save = '.stim2Impal\';
end

if ~isfield(opts, 'duration_sec')
    opts.duration_sec = 1.0;    % (sec)
end

if ~isfield(opts, 'nbits')
    opts.nbits = 16;    % (sec)
end

if ~isfield(opts, 'rms')
    opts.rms = 'win';    % (sec)
end


assert( isfield(opts, 'params'), '--> ERROR in [split_wav_to_Impale.m]: OPTS must contain <params>!' );
assert( isfield(opts.params, 'name'), '--> ERROR in [split_wav_to_Impale.m]: opts.params must contain <name>!' );
assert( isfield(opts.params, 'val'), '--> ERROR in [split_wav_to_Impale.m]: opts.params must contain <val>!' );

Fs = opts.Fs;

% Number of bits to save the WAV files
nbits = opts.nbits;     

% Create cell array to store filenames
file_list = {}; 

max_Y_len = cellfun(@(YY) length(YY), Y_list);
assert( 0 == sum(diff(max_Y_len)), '--> Error at [split_wav.m]: All entries must have the SAME LENGTH!');
max_Y_len = max_Y_len(1);


%% Maximum size of each segment (repetition)
% MAX_SEQ_SMP: The length of each chunk (for Impale)
max_seg_smp = fix(opts.duration_sec * Fs);   % (samples)
Nrep        = ceil(max_Y_len/max_seg_smp);
tot_len     = max_seg_smp*Nrep;
pad_zeros   = tot_len - max_Y_len;

rms_values = zeros(n_signals,2);
 
for row = 1:n_signals
    yij = Y_list{row};

    % The total stimulus; pad with zeros at the end
    ypad_ij = [yij; zeros(pad_zeros, size(yij,2))];

    % Left & right ears
    ypad_L = ypad_ij(:, 1);
    ypad_R = ypad_ij(:, 2);

    % Normalize; avoid data clippng when saving as WAV
    ypad_L = 0.999*ypad_L./max(abs(ypad_L));    
    ypad_R = 0.999*ypad_R./max(abs(ypad_R));    

    switch lower(opts.rms)
        case 'simple'    	% RMS (simple)
            rms_seg_L_db = db( rms(ypad_L) );
            rms_seg_R_db = db( rms(ypad_R) );
            
        case 'threshold'  	% RMS (global threshold)
            rms_thr = 0.1;  % avoid doing RMS over lower amplitude intervals
            rms_seg_L_db = db( rms(ypad_L(abs(ypad_L)>rms_thr,1)) );
            rms_seg_R_db = db( rms(ypad_R(abs(ypad_R)>rms_thr,1)) );
            
        case 'win'          % RMS (percentile)      
            rms_seg_L_db = rms_running_window(ypad_L, Fs);
            rms_seg_R_db = rms_running_window(ypad_R, Fs);
            
        otherwise
            error('--> [split_wav_to_Impale.m]: unrecognized <opts.rms> option!');
    end

    % Save the RMS values, if needed
    if 1 <= nargout
        rms_values(row,:) = [rms_seg_L_db, rms_seg_R_db];
    end
    
    % Create the filename with the parameter list
    param_list = arrayfun(@(X) num2str(X), opts.params.val(row,:), 'UniformOutput', false);
    param_list_Impale = strcat({'_-'}, opts.params.name, param_list);
    param_list_Impale = [param_list_Impale{:}];

    % Save the colmplete WAV file (i.e., without splitting it) for reference
    % =============================================================
    fn_R = [opts.path2save, opts.fn_generic, param_list_Impale, 'R'];
    fn_L = [opts.path2save, opts.fn_generic, param_list_Impale, 'L'];

    % Save the full & zero padded stimulus files
    epl.EPLwavwrite( ypad_L, Fs, nbits, fn_R, 'rms_dB', rms_seg_L_db );
    epl.EPLwavwrite( ypad_R, Fs, nbits, fn_L, 'rms_dB', rms_seg_R_db );

    % =============================================================

    % Repetitions: split the signal into repetirions
    if 1 < Nrep    % is it needed? if the starting segment is short (e.g., one 
                   % utterance) then there is no need for using the
                   % repetitions
        
        % Make sure that the path that holds all repetitions exists
        new_path = [opts.path2save, opts.fn_generic, param_list_Impale];
        if ~isdir(new_path)
            mkdir(new_path)
        end

        for tok = 1:Nrep   
            rep_idx = (tok-1)*max_seg_smp + (1:max_seg_smp);

            yrep_L = ypad_L(rep_idx);
            yrep_R = ypad_R(rep_idx);

            % Name stimulus -- RIGHT ear
            fn_tok_R = [new_path, filesep, opts.fn_generic, param_list_Impale,...
                '_-', 'Rep', num2str(tok), 'R'];

            % Name stimulus -- LEFT ear
            fn_tok_L = [new_path, filesep, opts.fn_generic, param_list_Impale,...
                '_-', 'Rep', num2str(tok), 'L'];

            % Save the stimulus file
            epl.EPLwavwrite( yrep_L, Fs, nbits, fn_tok_R, 'rms_dB', rms_seg_L_db );
            epl.EPLwavwrite( yrep_R, Fs, nbits, fn_tok_L, 'rms_dB', rms_seg_R_db );

            file_list = [file_list, {fn_tok_L}, {fn_tok_R}];
            %fprintf('--> [split_wav.m]: saved %s\n', fn_tok);
        end
    end

end















