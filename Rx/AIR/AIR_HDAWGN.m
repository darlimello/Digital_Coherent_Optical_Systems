function [AIR] = AIR_HDAWGN(SNRdB,ModFormat,AIREval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIRAWGN [AIR] = AIR_HDAWGN(SNRdB,ModFormat,AIREval)                     %
%                                                                         %
%   This function calculates approximations of the theoretical AIRs for   %
% schemes with HD-SW and HD-BW decoders under the circularly symmetric    %
% Gaussian channel. As symbol error probability (SEP) and bit error       %
% probability (BEP), we use theoretical pre-FEC error performance         %
% expressions. For BEP, we consider an error expression derived for Gray  %
% mapping. The supported modulation formats are 'QPSK' and '16QAM'.       %
%                                                                         %
% Input:                                                                  %
%   SNRdB     = SNR in dB;                                                %
%   ModFormat = Modulation format: 'QPSK' or '16QAM';                     %
%   AIREval   = Defines the AIR to the evaluated: 'HD-SW' or 'HD-BW';     % 
%                                                                         %
% Output:                                                                 %
%   AIR = Approximation of the theoretical AIR (in bit/symbol) for HD-SW  %
%         or HD-BW decoders (depending on the variable 'AIREval') under a %
%         circularly symmetric Gaussian channel.                          %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SNR:
    SNR = 10.^(SNRdB/10);

    % Defining parameters according to the modulation format:
    switch ModFormat
        case 'QPSK'
            % QPSK modulation
            if strcmp(AIREval,'HD-BW')
                BEP = (1/2)*erfc(sqrt(SNR/2));
            elseif strcmp(AIREval,'HD-SW')
                M = 4;
            end
        case '16QAM'
            % 16QAM modulation
            if strcmp(AIREval,'HD-BW')
                BEP = (3/8)*erfc(sqrt(SNR/10)) +...
                   (1/4)*erfc(3*sqrt(SNR/10)) - (1/8)*erfc(5*sqrt(SNR/10));
            elseif strcmp(AIREval,'HD-SW')
                M = 16;
            end
    end    
    if strcmp(AIREval,'HD-SW')    
        SEP = 2*(sqrt(M)-1)/sqrt(M)*erfc(sqrt(3*SNR/(2*(M-1))))...
              - (((sqrt(M)-1)/sqrt(M))*erfc(sqrt(3*SNR/(2*(M-1))))).^2;
    end
    
    % AIR estimate:
    switch AIREval
        case 'HD-SW'
            [AIR] = AIR_HDSW(SEP,ModFormat);
        case 'HD-BW'
            [AIR] = AIR_HDBW(BEP,ModFormat);
    end  
end