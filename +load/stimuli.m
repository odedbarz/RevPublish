function stim_st = stimuli(fs_resample, path2wav, fn_template, dist_list, revb_list)
%
% stim_st = load.stimuli(fs_resample, path2wav, fn_template,dist_list, revb_list)
%
% Description:
% Loads all stimuli from harddisk into the work space.   
%
%

import load.*

%% Set the input
if 1 > nargin
    fs_resample = 100e3;     %(kHz)
end

if 2 > nargin
    %path2wav = './.wav_files/Spch_(36)sec/';
    path2wav = load.path_to_data('wav_spch_36sec');
end

if 3 > nargin
    fn_template = 'BRIR_L_-Dist%d_-Rev%d.wav';
end


if 4 > nargin
    dist_list = [15 30];      % Dist. list
    revb_list = [2 8 10];     % Revb. list
end



%%
len_dist = length(dist_list);   % # Distance 
len_revb = length(revb_list);   % # Reverberation coefficient

% Create a table that will hold all information about the loaded stimuli
T = table([], [], 'VariableNames', {'Revb', 'Dist'});
for kk = 1:length(dist_list)
    T = [T; table(revb_list(:), dist_list(kk)*ones(size(revb_list(:))), 'VariableNames', {'Revb', 'Dist'})];
end
stim_st.table_meas = T;
stim_st.fn_template = fn_template;


%%
Y       = cell(len_revb, len_dist);
fn_mtx  = cell(len_revb, len_dist);

for jj = 1:len_dist
    for ii = 1:len_revb    
        fn_ij = sprintf(fn_template, dist_list(jj), revb_list(ii));
        full_fn = fullfile(path2wav, fn_ij);
        [Y{ii,jj}, info] = load.wav( full_fn, fs_resample );

        % Read the parameters from the file name
        fn_mtx{ii,jj} = fn_ij;
        Dist_ij = regexp(fn_mtx{ii,jj}, '(?<=Dist)\d+','once','match');
        assert(~isempty(Dist_ij), '--> ERROR @ [load_stimuli.m]: can''t find a Dist!!');
        Dist_ij = str2double(Dist_ij);
        
        Rev_ij = regexp(fn_mtx{ii,jj}, '(?<=Rev)\d+','once','match');
        assert(~isempty(Rev_ij), '--> ERROR @ [load_stimuli.m]: can''t find a Dist!!');
        Rev_ij = str2double(Rev_ij);
           
        % Set the labe of the file for later reference
        labels{ii,jj} = sprintf('%.1f m, %g%%', 1/10*Dist_ij, 10*Rev_ij);
        
        % Check that the duration of the signal is OK (legacy measurements of 36 sec)
        assert( 0==nnz( Y{ii,jj} > info.Duration ),...
            '--> ERROR: check out that you are using the right WAV files (the right length)!' );

    end
end
stim_st.labels = labels;


%% Set the output structure 
stim_st.Y     = Y;
stim_st.info  = info;
stim_st.fs    = stim_st.info.SampleRate; 

len_y = arrayfun(@(I) length(stim_st.Y{I}), 1:length(stim_st.Y(:)));
if 1 == length( unique(len_y) )
    stim_st.len = len_y(1);
    stim_st.t   = linspace(0, 1/stim_st.fs*stim_st.len, stim_st.len)'; 
else
    stim_st.len = len_y;
    stim_st.t   = arrayfun(@(X) linspace(0, 1/stim_st.fs, stim_st.len(X))',...
        len_y, 'UniformOutput', false );     
end
stim_st.units = 'sec';


































