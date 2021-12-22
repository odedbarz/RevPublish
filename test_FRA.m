%
% test_FRA.m
% 
% Description:
%   This script tests the FRA procedure. It finds the FRA rates, CF(s), and 
%   plot the results. 
% 
% Note: 
%   In order to run this script, you need to have an impale structure S on
%   the working memory!
%

clc
clear all

S = load('OBZ_C64_036-2-1-FRAL.mat');
% S = load('OBZ_C64_027-1-5-FRAL.mat');
% S = load('OBZ_C74_017-2-1-FRAL.mat');

viewer.plot_FRA(S);





