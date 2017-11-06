function RM = grouping(I, V, U, d)
%GROUPING combination of a subset of elementary matrices 
%
% Inputs:      I   --> Vector with elementary matrices indices
%              d   --> Eigenvalues
%              U   --> Eigenvectors
%              V   --> Factor Vectors
%
% Output:      RM  --> Resultant Matrix: sum of elementary matrices within a group
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The third step in SSA decomposition, performing the combination of a subset 
% of elementary matrices, that capture the desired structures. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Michel Leles, 07/05/2015

D=repmat(d', length(V), 1);

P=D.*V;

RM=U(:,I)*P(:,I)';



