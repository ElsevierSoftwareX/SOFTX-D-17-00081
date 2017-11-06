function [d,U,V] = svd_calc(A, method)
%SVD_CALC Compute Singular Value Decomposition
%
% Input:       A       --> Trajectory Matrix
%              method  --> SVD calculation method: 'bsc' = basic SSA and     
%                                                  'kor'= computational efficient
%
% Output:      d       --> eigenvalues
%              U       --> Eigenfilters
%              V       --> Factor Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second step in SSA decomposition, which factorizes the trajectory matrix 
% into a set of L elementary matrices: A_1, A_2, ..., A_L. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Michel Leles, 07/05/2015

if(strcmp(method, 'bsc'))
    [U,S,V] = svd(A);
elseif(strcmp(method, 'kor'))
    [U,S,V] = svdecon(A);
else % default
    [U,S,V] = svdecon(A);
end
d = diag(S);