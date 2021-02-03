function setup_environment(path_root)
%
%

if 0 == nargin
    path_root = pwd;
end

% Add the path of the StimViewerGUI GUI
addpath(path_root);

addpath(genpath(fullfile(path_root, 'CochlearModels')));

aux.figure_setup;


set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 