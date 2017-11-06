function y = diagonal_averaging(RM)
%DIAGONAL AVERAGING compute the diagonal averaging
%
% Input:      RM --> Resultant Matrix of Grouping Step
%
% Output:     y  --> Reconstructed time-series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The last step in SSA decomposition, transforming the resultant matrix from 
% grouping into the reconstructed time-series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Michel Leles, 07/05/2015

L=size(RM,1);
K=size(RM,2);
N=K+L-1;

y=zeros(N,1);
Lp=min(L,K);
Kp=max(L,K);

for k=0:Lp-2
    for m=1:k+1
        y(k+1)=y(k+1)+(1/(k+1))*RM(m,k-m+2);
    end
end

for k=Lp-1:Kp-1
    for m=1:Lp;
        y(k+1)=y(k+1)+(1/(Lp))*RM(m,k-m+2);
    end
end

for k=Kp:N
    for m=k-Kp+2:N-Kp+1
        y(k+1)=y(k+1)+(1/(N-k))*RM(m,k-m+2);
    end
end
%% Column vector
y=y';