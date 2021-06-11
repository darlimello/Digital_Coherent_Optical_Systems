function L = LLR(y,x,bmap,Sigma2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L] = LLR(y,x,bmap,TSigma)                                              %
%                                                                         %
%  This function calculated LLRs considering the exact expression for a   %
% circularly symmetric Gaussian channel, for the received symbols 'y',    %
% transmitted symbols 'x', and binary labels given by 'bmap'.             %
%                                                                         %
% Input:                                                                  %
%   y      = Input signal (one pol. orientation) normalized to unitary    %
%            power. 'y' must be a column vector;                          %
%   x      = Reference constellation normalized to unitary power. 'x' must%
%            be a column vector;                                          %
%   bmap   = Binary label of each symbol of the reference constellation   %
%            'x' (in corresponding order to the column vector 'x'). 'bmap'%
%            must be a matrix with 'm = log2(M)' columns and 'M' rows,    %
%            where 'M' is the modulation order;                           %
%   Sigma2 = Estimate of the variance of the circularly symmetric Gaussian%
%            channel. 'Sigma2' is the variance per dimension.             %
%                                                                         %
% Output:                                                                 %
%   L = Estimated LLRs (column vector);                                   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of bits per symbol:
    m = size(bmap,2);
    
    % LLR estimation:
    L = zeros(numel(y),m);    
    for k = 1:m        
        % Sets with bit 'b = {0,1}' at position 'k':
        xSet_b1 = x(bmap(:,k)==1) ; xSet_b0 = x(bmap(:,k)==0);
        
        % LLR estimation:
        num = zeros(numel(y),1) ; den = zeros(numel(y),1);        
        for i = 1:m^2/2          
            num = num + exp(-abs(y-xSet_b1(i)).^2/(2*Sigma2));
            den = den + exp(-abs(y-xSet_b0(i)).^2/(2*Sigma2));
        end                
        L(:,k) = log(num./den);         
    end
end