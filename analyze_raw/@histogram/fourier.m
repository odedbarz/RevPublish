function s = fourier(h, dim)
% HISTOGRAM/FOURIER - Compute Fourier spectrum of a histogram
% S = FOURIER(H) returns the complex spectrum (magnitude and phase) of H 
% as the spectrum object S.  
% For 2D histograms or neurograms, the spectrum is computed along each column of H, 
% except for raster histograms in which case the spectrum is computed along each row.
%
% FOURIER(H, DIM) returns the spectrum along dimension DIM.
%
    dmax = ndims(h.data);
    if nargin < 2,
        if strcmp(h.type, 'PESE') | ~isempty(strfind(h.type, 'Raster')),
            dim = dmax;
        elseif num_dim(h.type) == 2,
            dim = dmax-1;
        else
            dim = which_dim(h);
        end
    end
    nfft = size(h.data, dim);
    x = fft(full(h.data), nfft, dim);
    
    % remove negative frequencies
    for k = 1:dmax, subs{k} = ':'; end
    subs{dim} = 1:nfft/2+1;
    s = spectrum(x(subs{:})); 
    
    % create frequency vector in Hz
    f = 1000*[0:nfft/2]/(h.binwidth(dim-dmax+2)*nfft);
    
    % set variable names and values
    fname = 'Frequency (Hz)';
    s = set(s, 'VarName', fname, dim, fname, f);
    
    % specify second histogram dimension
    if num_dim(h.type) == 2,
        if dim == dmax, dim2 = dmax-1; else dim2 = dmax; end
        name = axis_label(h, dim2);
        s = set(s, 'VarName', name, dim2, name, get(h, 'BinTimes', dim2-dmax+2));
    end
    
    % set stimulus variables
    for k = 1:length(h.value),
        s = set(s, 'VarName', h.varname{k}, k, h.varname{k}, h.value{k});
    end
    
    