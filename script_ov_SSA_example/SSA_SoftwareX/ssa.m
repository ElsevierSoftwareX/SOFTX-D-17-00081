function y = ssa(series, L, I, method)
%SSA Wrapper function for SSA decomposition
%
% Inputs:      series --> Time-series under analysis
%              L      --> Embedding dimension
%              I      --> Vector with elementary matrices indices
%              method --> SVD calculation method: 'bsc' = basic SSA and     
%                                                 'kor'= computational efficient
%
% Output:      y      --> reconstrutcted time-series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SSA technique consist in four major steps. 
%
% Step 1: The Immersion, which maps the raw time-series into a Hankel matrix, 
% called trajectory matrix. Parameter L, the embedding dimension, is used to assemble the trajectory matrix. 
%
% Step 2: The Singular Value Decomposition (SVD), which factorizes the trajectory matrix 
% into a set of L elementary matrices: A_1, A_2, ..., A_L. 
%
% Step 3: In Grouping, occurs the combination of a subset of elementary matrices, that capture the desired structures. 
%
% Step 4: Diagonal Averaging transforms the resultant matrix from grouping into the reconstructed time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Michel Leles, 07/05/2015

% Step 1: Immersion  
A = buffer(series, L, L-1, 'nodelay');
%BUFFER Matlab function that buffer signal vector into matrix of data frames

% Step 2: SVD
[d,U,V] = svd_calc(A, method);

% Step 3: Groupíng
RS = grouping(I, V, U, d);

% Step 4: Diagonal Averaging
y =  diagonal_averaging(RS);