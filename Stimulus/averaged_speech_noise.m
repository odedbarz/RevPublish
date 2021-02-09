function [yavg, num_wav, Yfft] = averaged_speech_noise(Ywav, fignum)
% 
%    function yavg = shaped_speech_noise(Ywav, rand_opt)
%
% Ywav: a cell of stacked WAV files to average.
% 
% num_wav: number of files used for the average.
%
%
% Description:
% Average the spectrum of all the WAV files in Ywav. The phase is
% uniformly randomized between [0,2PI].
%
% 

    if 1 >= nargin, fignum = []; end;
        
    % Number of waves
    num_wav = length(Ywav);

    % The length of the output is the length of maximum signal
    len_avg = max(cellfun(@(X) length(X), Ywav));

    % Yavg = zeros(len_avg, Nwav);

    Yfft = cellfun(@(X) fft(X,len_avg), Ywav, 'UniformOutput', false);
    Yfft = cell2mat(Yfft);

    % Average spectrum, but without a phase (yet!)
    Y_avg_spectrum = mean(abs(Yfft),2);

    yavg = randomize_phase(Y_avg_spectrum);

    if ~isempty(fignum)
        figure(fignum);
        clf;
        plot(20*log10( abs(Yfft) ));
        title(sprintf('A batch of %d spectrum responses', num_wav));
        ylabel('Spectrum (dB)');
        xlabel('Samples');
        hold on
        plth = plot(20*log10( Y_avg_spectrum ), '--k');
        hold off
        legend(plth, 'Avg. Spectrum Response');
    end
    
end


function [yrphs, Yrphs] = randomize_phase(Y)
%
%
%
% Randomizes teh phase in the phase domain.
%

    isodd = size(Y,1)/2 ~= fix(size(Y,1)/2);

    if isodd    % Yv is ODD length
        dummy = Y(2:end);          % +1: don't take the DC
        Yv_half = dummy(1:end/2);
        Yv_dc   = Y(1);

        % Randomize the phase
        rnd_phs    = 2*pi*rand(size(Yv_half,1),1);     % the random phase
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
        rnd_phs    = 2*pi*rand(size(Yv_half,1),1);     % the random phase
        Yrphs_half = abs(Yv_half) .* exp(1j.*rnd_phs);

        % Stitch it all together
        Yrphs = [Yv_dc; Yrphs_half; Yv_middle; conj(Yrphs_half(end:-1:1))];
        yrphs = ifft( Yrphs ); %, 'symmetric' );
        assert(0 == sum(imag(yrphs)), '--> ERROR: you didn''t do a proper stitching!');

    end


end
