function [AIR] = AIR_HDBW(BEP,ModFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR_BWHD [AIR] = AIR_HDBW(BEP,ModFormat)                                %
%                                                                         %
%   This function estimates AIRs for CM schemes with hard-decision bit-   %
% wise (HD-BW) decoders. It is assumed that the decoding is based on      %
% the Hamming distance, so that we can interpret that the decoder sees 'm'%
% parallel binary symmetric channels with error probability 'BEP', where  %
% 'm = log2(M)' and 'M' is the modulation order. In fact, 'BEP' is the    %
% pre-FEC average error probability across the 'm' bit positions. 'BEP'   %
% can be estimated as the pre-FEC bit error rate (BER).                   %
%                                                                         %
% Input:                                                                  %
%   BEP       = Pre-FEC average bit error probability across the 'm' bit  %
%               positions. 'BEP' can be estimated as the pre-FEC BER (per %
%               pol. orientation);                                        %
%   ModFormat = Modulation format: 'QPSK' or '16QAM';                     %
%                                                                         %
% Output:                                                                 %
%   AIR = AIR estimate for a CM scheme with HD-BW decoding;               %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Defining parameters according to the modulation format:
    switch ModFormat  
        case 'QPSK'            
            M = 4;
        case '16QAM'
            M = 16;
    end
    % Binary entropy calculated with 'BEP' and AIR estimate:
    hb  = -BEP.*log2(BEP) - (1-BEP).*log2(1-BEP);
    AIR = log2(M)*(1 - hb);        
end