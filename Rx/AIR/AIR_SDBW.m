function [AIR] = AIR_SDBW(x,b,y,ModFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR_SDBW [AIR] = AIR_SDBW(x,b,y,ModFormat)                              %
%                                                                         %
%  This function estimates AIRs for CM schemes with soft-decision bit-wise%
% (SD-BW) decoders. It is assumed a circularly symmetric Gaussian channel:%
%                               y = h*x + n                               %  
% where 'y' is the received symbol, 'x' is the transmitted symbol, 'h' is %
% a scalar and 'n' is a complex Gaussian random variable with total       %
% variance TSigma2 = 2*sigma^2. The AIRs are estimated using LLRs. LLRs   %
% are calculated using the exact expression for a circularly symmetric    %
% Gaussian channel. Gray mapping is assumed. All input sequences must     %
% be synchronized.                                                        %
%                                                                         %
% Input:                                                                  %
%   x         = Sequence of symbols transmitted in one pol. orientation   %
%               (column vector) normalized to unitary power;              %
%   b         = Sequence of bits transmitted in one pol. orientation      %
%               (column vector) normalized to unitary power;              %
%   y         = Sequence of symbols received in one pol. orientation      %
%               (column vector);                                          %
%   ModFormat = Modulation format;                                        %
%                                                                         %
% Output:                                                                 %
%   AIR = AIR estimate for a CM scheme with BW-SD decoding;               %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Defining parameters according to the modulation format:
    switch ModFormat  
        case 'QPSK'
            % Modulation order and number of bits per symbol:
            M = 4 ; m = log2(M);
            % Bit-mapping:
            bmap = [0 1; 0 0; 1 1; 1 0];
            % Constellation symbols:
            s = [-1+1i; 1+1i; -1-1i; 1-1i].' ; s = s/sqrt(2);
        case '16QAM'
            % Modulation order and number of bits per symbol:
            M = 16 ; m = log2(M);
            % Bit-mapping:
            bmap = [0 0 0 0; 0 1 0 0; 0 1 0 1; 0 0 0 1; 1 0 0 0;...
                    1 1 0 0; 1 1 0 1; 1 0 0 1; 1 0 1 0; 1 1 1 0;...
                    1 1 1 1; 1 0 1 1; 0 0 1 0; 0 1 1 0; 0 1 1 1; 0 0 1 1];
            % Constellation symbols:          
            s = [3+3i; 1+3i; -1+3i; -3+3i; 3+1i; 1+1i; -1+1i; -3+1i;...
                 3-1i; 1-1i; -1-1i; -3-1i; 3-3i; 1-3i; -1-3i; -3-3i].';
            s = s/sqrt(10);
    end

    % Estimating the variance of the complex Gaussian distr.:
    h = x'*y/norm(x)^2 ; Sigma2 = mean(abs(y-h*x).^2)/2;
    
    % Calculating LLRs:
    L = LLR(y,h*s,bmap,Sigma2);
    
    % Reshaping the transmitted bits:
    b = reshape(b',m,numel(b)/m)';

    % Computing the AIR (performing minimization over 'a'):
    [~,AIRaux] = fminbnd(@(a) sum(mean(log2(1 + exp(a*(-1).^b.*L)))),0,2);
    AIR        = m - AIRaux;    
end