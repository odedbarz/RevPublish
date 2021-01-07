function Xx = compute_DFT_via_Goertzel(xin,f,Fs)
%
% This function was taken from MAT:AB's signal processing toolbox, 
%   C:\Program Files\MATLAB\R2017a\toolbox\signal\signal\private\
% 
% I use it specifically for computinf DFT for nonuniform spaced
% frequencies, e.g., exponential scale frequencies.
%
%
% Use Goertzel to compute raw DFT

f = mod(f,Fs);    % 0 <= f < = Fs                                                    
xm = size(xin,1); % NFFT

% wavenumber in cycles/period used by the Goertzel function
% (see equation 11.1 pg. 755 of [2]) 
k = f/Fs*xm;

% call Goertzel for each channel
Xx = zeros(size(k,1),size(xin,2));
for i=1:size(xin,2)
    % goertzelmex is a MEX function taken from MATLAB's signal processing
    % toolbox    
    %Xx(:,i) = spec.goertzelmex(double(xin(:,i)),k);
    
    computer_name = getenv('computername');
    if strcmpi('OdedAlienware', computer_name)
        Xx(:,i) = spec.goertzelmex1(double(xin(:,i)),k);

    elseif strcmpi('ODEDSDELL', computer_name)
        % Oded's external HD
        Xx(:,i) = spec.goertzelmex2(double(xin(:,i)),k);

    else        
        % Apollo drive
        error('[compute_DFT_via_Goertzel]: You need to have the right MEX file!');
    end

    
    
end

