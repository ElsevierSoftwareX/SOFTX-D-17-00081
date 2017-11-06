function y = segmentation(x, L, I, Z, q, method)
%SEGMENTATION compute the ov-SSA algorithm
% Inputs:       x      --> Original time-series
%               L      --> Embedding dimension
%               I      --> Vector with elementary matrices indices
%               Z      --> Local segment lenght
%               q      --> The number of samples that are reconstructed locally
%               method --> SVD calculation method: 'bsc' = basic SSA and     
%                                                  'kor'= computational efficient
%
% Outputs:      y      --> Reconstructed SSA time-series
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider a time series segment of length Z. 
% All samples within this segment are used to compute the SSA locally. 
% However, since the SSA method suffers from boundary effects, 
% the extreme points in the left and right edges are discarded. 
% The quantity of discarded samples is given by L_B = (Z-q)/2. 
% Only an inner subset of samples, q, is considered meaningful to represent the local time-series Z. 
% The final reconstruction is given by the concatenation of the inner segmentes q, which do not overlap. 
% The extreme edges of the original time-series need special attention.
% In the first run only the last L_B points are discarded.
% On the other hand, in the last run, the first L_B points are discarded.
% This approach is a modification of  the overlap-save method, 
% a classic tool to calculate FFT (Fast Fourier Transform) convolution of 
% infinite (and finite) duration time series. 
% This adaptation was necessary because the standard SSA algorithm suffers 
% from boundary effects on both sides. In the overlap-save method only the initial 
% points must be discarded, because boundary effects occur only at the filtering initialization. 
% For a complete discussion about the method, see Leles et al (2017) A New
% Algorithm in Singular Spectrum Analysis Framework: the Overlap-SSA (ov-SSA). 
% SoftwareX. In Press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Michel Leles, 07/05/2015

%% Parameters of algorithm
N=length(x); % length of time-series
y=zeros(1,N);
P=floor((N-Z)/q)+1; % the number of iterations
L_B=(Z-q)/2; % number of points discarded at each iteration
%% L_B must be an integer, otherwise a Matlab error is launched
if(mod(L_B,1)~=0)
    error('The number of points discarded at each iteration must be an integer')
end

%% First run
p=1;
series = x(p:p+Z-1);
y_aux=ssa(series, L, I, method);
y(1:Z-L_B)=y_aux(1:Z-L_B); %
%% Loop

for p=2:P-1
    if(p==1)  % first run
        
    else % concatenation
        rho=(p-1)*q+1;
        series = x(rho:rho+Z);
        y_aux=ssa(series, L, I, method);
        y(rho+L_B:rho+L_B+q-1)=y_aux(L_B+1:L_B+q)';
    end
end
%% Last iteration
p=p+1;
rho=(p-1)*q+1;
series = x(rho:end);
y_aux=ssa(series, L, I, method);
y(rho+L_B:end)=y_aux(L_B+1:end)';