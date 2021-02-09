function [text_h, ax, tbl_phonemes, tbl_metadata] = plot_phonemes(sp, varargin)
%
%   function [text_h, ax, tbl_phonemes, tbl_metadata] = plot_phonemes(sp, varargin)
%





%% Parse the input
p = inputParser;

addRequired(p, 'sp', @isnumeric);

addOptional(p, 'duration_sec', 36, @isnumeric);      % (sec)
addOptional(p, 'tbl_metadata', [], @istable);
addOptional(p, 'ax', [], @ishandle);
addOptional(p, 'fontsize', 22, @isnumeric);
addOptional(p, 'verbose', 0, @isnumeric);
addOptional(p, 'fignum', 0, @isnumeric);

parse(p, sp, varargin{:});

pars = p.Results;




%% Load the the META-DATA of the stimuli (TIMIT) files
if isempty(pars.tbl_metadata)
    %duration_sec = 1e-3*data.stim_st.duration_ms;
    fn_path_meta = load.path_to_data('Stimulus');
    fn_file_meta = sprintf('metadata_(%d)_wav_(30-Jun-2020)', pars.duration_sec);
    dummy        = load( fullfile( fn_path_meta, fn_file_meta ) );
    tbl_metadata = dummy.tbl_metadata;
else
    tbl_metadata = pars.tbl_metadata;
end

if isempty(pars.ax) || ~ishandle(pars.ax)
    ax = gca; 
else
    ax = pars.ax;
end




%%
% New table of phonemes of a selected speaker
tbl_phonemes = [table(tbl_metadata.phn{sp}, 'VariableNames', {'phn'}),...
    table(tbl_metadata.phn_start{sp} - tbl_metadata.t0_timit(sp), 'VariableNames', {'phn_start'}),...
    table(tbl_metadata.phn_end{sp} - tbl_metadata.t0_timit(sp), 'VariableNames', {'phn_end'}),...
];

tbl_phonemes = [tbl_phonemes, ...
    table(tbl_phonemes.phn_start * 1/tbl_metadata.fs(sp), 'VariableNames', {'t0'}),...
    table(tbl_phonemes.phn_end * 1/tbl_metadata.fs(sp), 'VariableNames', {'t1'})...
];


% PLOT
set(ax, 'XTickLabels', '');
set(ax, 'YTickLabels', '');
ylim([0, 1]);

aux.vline(tbl_phonemes.t0(1), 'ax', ax, 'Color', 'k');
aux.vline(tbl_phonemes.t1, 'ax', ax, 'Color', 'k');
linkaxes(ax, 'x');

d_times = tbl_phonemes.t1 - tbl_phonemes.t0;
t_phn = tbl_phonemes.t0 + 0.5*d_times;
phn_for_text = aux.mName2latex(tbl_phonemes.phn);
text_h = text(ax, t_phn, mean(ylim)*ones(length(t_phn),1), phn_for_text,...
    'HorizontalAlignment', 'center', 'FontSize', pars.fontsize);













