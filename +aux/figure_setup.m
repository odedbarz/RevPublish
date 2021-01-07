function opts = FigSetup
%
% Set default figure options
% 
% Created by Oded Barzelay, 02/22/2018
%
% Notes:
%   Get factory settings in MATLAB:
%     diary('c:\temp\junk.txt'); % Save long output of next command
%     get(groot,'factory'); % List all factory settings
%     edit('c:\temp\junk.txt'); % Look for font settings
%


% opts.fig_flag   = 1;                                        % plot?
%opts.fignum     = 99;

opts.ms         = 20;                                       % markersize;
opts.dotraster.ms=10;
opts.marker     = '.';
opts.lw         = 2;                                        % linewidth;
opts.ls         = '-';                                      % linestyle       
opts.mkface     = 'flat';                                   % MarkerFaceColor
opts.dotcolors  = {aux.rpalette('scarlet'), aux.rpalette('blue')};  % (2x2) dot colors 

% Patch
opts.patch.color = 'red';
opts.patch.faceAlpha = 0.2;
opts.patch.lineStyle = 'None';

opts.xscale     = 'log';    % {'log', 'linear'}

% fontsize
opts.fontsize.text  = 14;                                      
opts.fontsize.axes  = 14;                                       
opts.fontsize.title = 18;                                       

opts.fontsize.titleMultiplier = 1.4;                                       % title size


%% Set MATLAB's figure defaults:
set(0, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter','latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');

set(groot, 'defaultLineLineWidth', opts.lw );
set(groot, 'defaultLineMarkerSize', opts.ms);

% Text Properties:
set(groot, 'DefaultTextFontSize', opts.fontsize.text);

% Axes properties:
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', opts.fontsize.axes, ...
      'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
      'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
      'DefaultAxesTitleFontWeight', 'normal', ...
      'DefaultAxesTitleFontSizeMultiplier', opts.fontsize.titleMultiplier ) ;


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  