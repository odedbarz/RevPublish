function yrphs = shaped_speech_noise(y, rand_opt)
% 
%    function yrphs = shaped_speech_noise(y, rand_opt)
%
% rand_opt == 1: random phase
% rand_opt == 2: shuffled phase
%
%
% Description:
% Speech spectrum with random phase + AM with the envelope, that is, Keep the speech
% spectrum but randomize the phase.

if 1 >= nargin, rand_opt = 1; end

Y = fft(y);

isodd = size(y,1)/2 ~= fix(size(y,1)/2);

if isodd    % Yv is ODD length
    dummy = Y(2:end);          % +1: don't take the DC
    Yv_half = dummy(1:end/2);
    Yv_dc   = Y(1);
    
    % Randomize the phase
    switch rand_opt
        case 1
            rnd_phs = 2*pi*rand(size(Yv_half,1),1);     % the random phase
        case 2
            shuffled_phases = randperm(size(Yv_half,1));
            rnd_phs = 2*pi*angle(Yv_half(shuffled_phases));     % the random phase
        otherwise
            error('--> Error in [shaped_speech_noise.m]: Unrecognized <rand_opt> option!');
    end
    Yrphs_half = abs(Yv_half) .* exp(1j.*rnd_phs);

    % Stitch it all together
    Yrphs = [Yv_dc; Yrphs_half; conj(Yrphs_half(end:-1:1))];
    yrphs = ifft( Yrphs ); %, 'symmetric' );
    assert(0 == sum(imag(yrphs)), '--> ERROR: you didn''t do a proper stitching!');

else    % Yv is EVEN length
    Yv_half   = Y(2:end/2);          % +1: don't take the DC
    Yv_dc     = Y(1);
    Yv_middle = Y(end/2+1);
    assert(isreal(Yv_middle), '--> ERROR: for an EVEN length signal, the center DFT coefficient is always real!');
    
    % Randomize the phase
    %rnd_phs    = 2*pi*rand(size(Yv_half,1),1);     % the random phase
    switch rand_opt
        case 1
            rnd_phs = 2*pi*rand(size(Yv_half,1),1);     % the random phase
        case 2
            shuffled_phases = randperm(size(Yv_half,1));
            rnd_phs = 2*pi*angle(Yv_half(shuffled_phases));     % the random phase
        otherwise
            error('--> Error in [shaped_speech_noise.m]: Unrecognized <rand_opt> option!');
    end    
    Yrphs_half = abs(Yv_half) .* exp(1j.*rnd_phs);

    % Stitch it all together
    Yrphs = [Yv_dc; Yrphs_half; Yv_middle; conj(Yrphs_half(end:-1:1))];
    yrphs = ifft( Yrphs ); %, 'symmetric' );
    assert(0 == sum(imag(yrphs)), '--> ERROR: you didn''t do a proper stitching!');
    
end

