% main_anova.m
% 
% Description:
% Perform ANOVA analysis for the DRRs vs SUs & MUAs; which one is "better"?
%


close all;
clc;
clear all;


fignum = 10;
verbose = 1;

setup_environment('../');




%% Testing Two-Ways ANOVA
% ref: https://www.mathworks.com/help/stats/two-way-anova.html
%

drr = get_DRR_list_and_indices;

% Lets create a dummy data base:
n_rep = 30;     % number of repetitions
n_drr = drr.n_drr;
A = randn(n_rep, n_drr);   

bias = 1;
for k = 1:n_drr
    A(:,k) = A(:,k) + k*bias;    
end

figure(1); clf;
plot_dotbox(A);

[p, tbl, stats] = anova2(A, n_rep, 0);
tbl = cell2table(tbl);
tbl

figure(2); clf;
c = multcompare(stats);
c


%%