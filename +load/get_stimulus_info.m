function fn = get_stimulus_info(varargin)
%
% function fn = get_stimulus_info([S] or [suffix], [duration_sec])
%
% Description:
%   Get the path & file templates for a given stimulus.
%

fn.path = '';
fn.template = '';

if 1 == length(varargin)
    S = varargin{1};
    suffix       = S.measParam.Suffix;
    duration_sec = S.stimChans{1}.Source.numTokens;    
elseif 2 == length(varargin)
    suffix       = varargin{1};
    duration_sec = varargin{2};    
else
    error('--> ERROR in [get_stimulus_info] wrong number of inputs!!');
end

% Select the right stimulus type
switch lower(suffix)
    case 'spch'
        % I switched bectween two stimuli; that's why I need this IF
        % statement
        if  40 == duration_sec 
            fn.path = load.path_to_data('wav_spch_40sec');

        elseif 36 == duration_sec  
            fn.path = load.path_to_data('wav_spch_36sec');

        else
            error('--> ERROR in [main_analuze_measurements.m]: Unrecognized <S.dacInt.numTotalTrials>!');

        end

        fn.template = 'BRIR_L_-Dist%d_-Rev%d.wav';

    
    case 'ns_spch_fc4khz'
        % I switched bectween two stimuli; that's why I need this IF
        % statement
        if  40 == duration_sec 
            fn.path = load.path_to_data('wav_nsspch_fc(1000hz)_40sec');            
        else
            error('--> ERROR in [get_stimulus_info.m]: Invalid stimulus!');            
        end

        fn.template = 'ns_Spch_fc(4000Hz)_L_-Dist%d_-Rev%d.wav';

    
    case 'ns_spch_fc1khz'
        % I switched bectween two stimuli; that's why I need this IF
        % statement
        if  40 == duration_sec  
            fn.path = load.path_to_data('wav_nsspch_fc(4000hz)_40sec');            
        else
            error('--> ERROR in [get_stimulus_info.m]: Invalid stimulus!');            
        end

        fn.template = 'ns_Spch_fc(1000Hz)_L_-Dist%d_-Rev%d.wav';

    
    otherwise
        error('--> ERROR: unrecognized stimulus type!!!');

end

fn.isdir = isdir(fn.path);









