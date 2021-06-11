function [BER] = ErrorCounting_DiffEnc(Rx,TxBits,ModFormat,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ErrorCounting [BER] = ErrorCounting_DiffEnc(Rx,Tx,ModFormat,...         %
%                                 MapperMode,NPol)                        %
%                                                                         %
%   This function evaluates the BER using the transmitted bit sequence    %
% 'TxBits' and the received symbol sequence 'Rx', considering differential%
% encoding. The received sequence of bits is generated after applying     %
% differential decoding on the sequence 'Rx'. The sequence 'Rx' and the   %
% sequence of transmitted symbols associated with the bit sequence        %
% 'TxBits' (i.e. 'Tx') must be synchronized prior to error evaluation. The%
% sequence 'Rx' must be normalized to unitary power.                      %
%                                                                         %
% Input:                                                                  %
%   Rx        = Sequence of received symbols normalized to unitary power*.%
%   TxBits  = Sequence of transmitted bits associated with sequence 'Tx'*;%
%   ModFormat = Modulation format: 'QPSK' or '16-QAM';                    %
%   NPol      = Number of pol. orientations used;                         %
%  *Note: For transmission in single pol. orientation, 'TxBits' and 'Rx'  %
%         must be column vectors. For transmission with pol. multiplexing,%
%         'TxBits' and 'Rx' must be matrices with two column vectors,     %
%         where each column vector corresponds to the signal of one pol.  %
%         orientation.                                                    %
%                                                                         %
% Output:                                                                 %
%   BER = Struct with the evaluated BER per pol. orientation. If NPol == 1%
%         the struct BER has only the field BER.V. If NPol == 2, the      %
%         struct BER has the filed BER.V and BER.H;                       %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Decision of the received symbols V polarization:
    DRxBitsV = Decision_DiffEnc(Rx(:,1),ModFormat);

    % Counting the bit error rate:
    BER.V = mean(DRxBitsV ~= TxBits(:,1));

    % If pol. multiplexing is used:
    if NPol == 2                 
        % Decision of the received symbols in H polarization:
        DRxBitsH = Decision_DiffEnc(Rx(:,2),ModFormat);
       
        % Counting the bit error rate:
        BER.H = mean(DRxBitsH ~= TxBits(:,2));
    end    
end