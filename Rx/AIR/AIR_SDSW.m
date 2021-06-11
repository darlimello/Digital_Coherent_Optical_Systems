function [AIR] = AIR_SDSW(x,y,ModFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR_SDSW [AIR] = AIR_SDSW(x,y,ModFormat)                                %
%                                                                         %
%   This function estimates AIRs for CM schemes with soft-decision        %
% symbol-wise (SD-SW) decoders. It is assumed a circularly symmetric      %
% Gaussian channel                                                        %
%                              y = h*x + n                                %  
% where 'y' is the received symbol, 'x' is the transmitted symbol, 'h' is %
% a scalar and 'n' is a complex Gaussian random variable with total       %
% variance TSigma = 2*sigma^2. All input sequences must be synchronized.  %
%                                                                         %
% Input:                                                                  %
%   x         = Sequence of symbols transmitted in one pol. orientation   %
%               (column vector) normalized to unitary power;              %
%   y         = Sequence of symbols received in one pol. orientation      %
%               (column vector) normalized to unitary power;              %
%   ModFormat = Modulation format: 'QPSK' or '16QAM';                     %
%                                                                         %
% Output:                                                                 %
%   AIR = AIR estimate for a CM scheme with SW-SD decoding;               %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Defining parameters according to the modulation format:
    switch ModFormat  
        case 'QPSK'
            % Modulation order and constellation symbols:
            M = 4 ; s = [-1+1i; 1+1i; -1-1i; 1-1i] ; s = s/sqrt(2);
        case '16QAM'
            % Modulation order and constellation symbols:
            M = 16 ; s = [3+3i; 1+3i; -1+3i; -3+3i; 3+1i; 1+1i; -1+1i; ...
                -3+1i; 3-1i; 1-1i; -1-1i; -3-1i; 3-3i; 1-3i; -1-3i; -3-3i];
            s = s/sqrt(10);
    end

    % Estimating the variance of the complex Gaussian distr.:
    h = x'*y/norm(x)^2 ; Sigma2 = mean(abs(y-h*x).^2)/2;
    
    % Estimation of the conditional pdf qY|X(y|x):
    qY_X = (1/(2*pi*Sigma2)*exp(-abs(y-h*x).^2/(2*Sigma2)));
    
    % Estimation of the conditional pdf qY|X(y|s), where 's' represents 
    % constellation symbols:
    qY_S = sum((1/(2*pi*Sigma2)*exp(-abs(repmat(y,1,numel(s))-h.*...
        repmat(s.',size(y,1),1)).^2/(2*Sigma2))),2);
    
    % Estimating the AIR:
    AIR = sum(log2(qY_X./(1/M*qY_S)))/numel(y);
end