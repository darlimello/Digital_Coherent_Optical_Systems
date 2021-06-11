function [AIR] = AIR_HDSW(SEP,ModFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR_SWHD [AIR] = AIR_HDSW(SEP,ModFormat)                                %
%                                                                         %
%   This function estimates AIRs for CM schemes with hard-decision symbol-%
% wise (HD-SW) decoders. It is assumed that the decoding is based on      %
% the Hamming distance, so that we can interpret that the decoder sees a  %
% 'M'-ary symmetric channel with error probability 'SEP', where 'M' is the%
% modulation order. In fact, 'SEP' can be estimated as the pre-FEC symbol %
% error rate (SER).                                                       %
%                                                                         %
% Input:                                                                  %
%   SEP       = Pre-FEC symbol error probability. 'SEP' can be estimated  %
%               as the pre-FEC SER (per pol. orientation);                %
%   ModFormat = Modulation format;                                        %
%                                                                         %
% Output:                                                                 %
%   AIR = AIR estimate for a CM scheme with HD-SW decoders;               %
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
    % Binary entropy calculated with 'SEP' and AIR estimate:
    hb  = -SEP.*log2(SEP) - (1-SEP).*log2(1-SEP);
    AIR = log2(M) - hb - SEP*log2(M-1);        
end