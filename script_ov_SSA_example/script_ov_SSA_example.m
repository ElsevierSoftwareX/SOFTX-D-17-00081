% Run this example to see some of the method funcionalities
clear all
close all
clc
%%
addpath([pwd,'/SSA_SoftwareX']);
%% Real world time series chosen 
% Monthly accidental deaths in the USA, between 1973 and 1978. 
x=load('USA_Death.dat')'; % column vector
% This time series was used in different SSA studies and is freely available for download at: 
% https://datamarket.com/data/set/22p0/accidental-deaths-in-usa-monthly-1973-1978
N=length(x);
%% SSA parameters - following Hassani (2007). Singular Spectrum Analysis:
%          Methodology and Comparison. Journal of Data Science, 5, 239-257. 
L=24;
I=1:13;
%% SSA calculation
y_ssa=ssa(x, L, I, 'kor');
%% ov-SSA parameters
Z=60;
q=2; 
% Attention: the number of points discarded at each iteration (L_B=(Z-q)/2), 
% must be an integer, otherwise a Matlab error is launched.
%% ov-SSA calculation
y_ov_SSA = segmentation(x, L, I, Z, q, 'kor');
%% Mean Absolute Error - a MATLAB anonymous function 
MAE = @(x,y) sum(abs(x-y))/length(x);
%% Reconstructed Time series (MAE) error 
MAE_SSA=MAE(x, y_ssa)
MAE_0v_SSA=MAE(x, y_ov_SSA)