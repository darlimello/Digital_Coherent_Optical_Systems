function [BER,SER] = ErrorCounting(Rx,Tx,TxBits,ModFormat,NPol,FPhaseTest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ErrorCounting [BER,SER] = ErrorCounting(Rx,Tx,TxBits,ModFormat,NPol,... %
%                                 FPhaseTest)                             %
%                                                                         %
%   This function evaluates the bit error rate (BER) and the symbol error %
% rate (SER) using the transmitted 'Tx' and received 'Rx' sequences. For  %
% BER evaluation, the received bit sequence is compared to the tranmitted %
% bit sequence 'TxBits'. The received bit sequence is generated after hard%
% decide the received symbols and get their binary labels. For SER        %
% evaluation, symbols of sequence 'Rx' are hard-decided and compared to   %
% the symbols of sequence 'Tx'. Both sequences 'Rx' and 'Tx' must be      %
% normalized to unitary power and synchronized prior to error evaluation. %
% This function considers Gray mapping, with the labels of the function   %
% 'SymbolGeneration'.                                                     %
%                                                                         %
% Input:                                                                  %
%  Rx         = Sequence of received symbols normalized to unitary power*.%
%  Tx         = Sequence of transmitted symbols normalized to unitary     %
%               power*.                                                   %
%  TxBits     = Sequence of transmitted bits associated with sequence     %
%               'Tx'*;                                                    %
%  ModFormat  = Modulation format: 'QPSK' or '16-QAM';                    %
%  NPol       = Number of pol. orientations used;                         %
%  FPhaseTest = Flag to enable ('true') or disable ('false') the          %
%               evaluation of BER/SER for 4 different rotations (0 pi/2   %
%               pi 3*pi/2) applied to the entire sequence 'Rx'. This can  %
%               can be helpful for scenarios where 'Rx' may be a version  %
%               of 'Tx' simply rotated by a multiple of pi/2 (due to the  %
%               phase ambiguity of the constellations), e.g., in simul.   %
%               with phase noise and phase recovery;                      %
% *Note: For transmission in single pol. orientation, sequences 'TxBits', %
%        'Tx', and 'Rx' must be column vectors. For transmission with pol.%
%        multiplexing, sequences 'TxBits', 'Tx', and 'Rx' must be matrices%
%        with two column vectors, where each column vector corresponds to %
%        the sequence transmitted in one pol. orientation.                %
%                                                                         %
% Output:                                                                 %
%   BER = Struct with the evaluated BER per pol. orientation. If NPol == 1%
%         the struct BER has only the field BER.V. If NPol == 2, the      %
%         struct BER has the filed BER.V and BER.H;                       %
%   SER = Struct with the evaluated SER per pol. orientation. If NPol == 1%
%         the struct SER has only the field SER.V. If NPol == 2, the      %
%         struct SER has the filed SER.V and SER.H;                       %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Test phase:
    if FPhaseTest
        RotAngle = [0 pi/2 pi 3*pi/2];
    else
        RotAngle = 0;
    end

    % Initializing an auxiliary vector:
    BERTempV = NaN(1,length(RotAngle)); SERTempV = NaN(1,length(RotAngle));

    for ii = 1:length(RotAngle)
        % Test phase:
        if FPhaseTest
            RxRot = Rx(:,1)*exp(1i*RotAngle(ii));
        else
            RxRot = Rx(:,1);
        end
        
        % Decision of the received symbols V polarization:
        DRxBitsV = Decision(RxRot,ModFormat,true);
        DRxSymbV = Decision(RxRot,ModFormat,false);

        % Counting the bit error rate:
        BERTempV(ii) = mean(DRxBitsV ~= TxBits(:,1));
        SERTempV(ii) = mean(DRxSymbV ~= Tx(:,1));
    end

    % Minumum value of BER:
    BER.V = min(BERTempV) ; SER.V = min(SERTempV);

    % If pol. multiplexing is used:
    if NPol == 2
        % Initializing an auxiliary vector:
        BERTempH=NaN(1,length(RotAngle)); SERTempH=NaN(1,length(RotAngle));
        
        for ii = 1:length(RotAngle)
          % Test phase:
          if FPhaseTest
               RxRot = Rx(:,2)*exp(1i*RotAngle(ii));
          else
               RxRot = Rx(:,2);
          end
          
          % Decision of the received symbols in H polarization:
          DRxBitsH = Decision(RxRot,ModFormat,true);
          DRxSymbH = Decision(RxRot,ModFormat,false);

          % Counting the bit error rate:
          BERTempH(ii) = mean(DRxBitsH ~= TxBits(:,2));
          SERTempH(ii) = mean(DRxSymbH ~= Tx(:,2));
        end

        % Minumum value of BER:
        BER.H = min(BERTempH) ; SER.H = min(SERTempH); 
    end    
end