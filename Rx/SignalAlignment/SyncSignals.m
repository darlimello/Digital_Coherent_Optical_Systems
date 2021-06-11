function [TxBits,Tx,Rx] = SyncSignals(TxBits,Tx,Rx,ModFormat,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SyncSignals [TxBits,Tx,Rx] = SyncSignals(TxBits,Tx,Rx,ModFormat,NPol)   %
%                                                                         %
%   This function synchronizes the sequence of transmitted symbols 'Tx' to%
% the sequence of received symbols 'Rx'. If needed, symbols of the        %
% sequence 'Tx' are discarded so that, at the output, both sequences have %
% the same length. Note: Symbols of sequence 'Rx' may also be discarded to%
% guarantee 'Tx' and 'Rx' are the same length at the output. The sequence %
% of transmitted bits associated with the 'synchronized' sequence 'Tx' is %
% also generated.                                                         %
%                                                                         %
% Input:                                                                  %
%   TxBits = Sequence of transmitted bits associated with sequence 'Tx'*; %
%   Tx     = Sequence of transmitted symbols normalized to unitary power*;%
%   Rx     = Sequence of received symbols normalized to unitary power*;   %
%   NPol   = Number of pol. orientations used;                            % 
%   *Note  = For transmission in single pol. orientation, 'TxBits', 'Tx', %
%            and 'Rx' must be column vectors. For transmission with pol.  %
%            multiplexing, 'TxBits', 'Tx', and 'Rx' must be matrices with %
%            two column vectors, where each column vector corresponds to  %
%            the signal of one pol. orientation;                          %
%                                                                         %
% Output:                                                                 %
%   Tx,Rx = Synchronized sequences. Sequences 'Tx' and 'Rx' are the same  %
%           length.                                                       %
%   TxBits = Sequence of transmitted bits associated with sequence 'Tx';  % 
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Paremeters related to the modulation format:
    switch ModFormat
        case 'QPSK'
            m = 2 ; AbsCorr = false;
        case '16QAM'
            m = 4 ; AbsCorr = true;
    end

    % Auxiliary variables:
    if NPol == 2
        Aux = 4;
    else
        Aux = 1;
    end
    Ratio = zeros(Aux,1) ; Pos   = zeros(Aux,1);

    % Cross-corr. between Tx 1 and Rx 1:
    [Pos(1),Ratio(1)] = Corr(Tx(:,1),Rx(:,1),AbsCorr);
            
    % If pol. multiplexing is used:
    if NPol == 2
       % Cross-corr. between Tx 2 and Rx 1:
       [Pos(2),Ratio(2)] = Corr(Tx(:,2),Rx(:,1),AbsCorr);
         
       % Cross-corr. between Tx 2 and Rx 2:
       [Pos(3),Ratio(3)] = Corr(Tx(:,2),Rx(:,2),AbsCorr);
        
       % Cross-corr. between Tx 1 and Rx 2:
       [Pos(4),Ratio(4)] = Corr(Tx(:,1),Rx(:,2),AbsCorr);
    end 

    % Checking if the pol. orientations are switched:
    if NPol == 2
        AuxFlip = false;
        if Ratio(2) >= Ratio(1) && Ratio(4) >= Ratio(3)          
            AuxFlip = true ; Aux     = [1 3];
        elseif Ratio(2) >= Ratio(1) && Ratio(3) >= Ratio(4)
            % Signals of pol. V and H may be switched:
            if Ratio(2) >= Ratio(3)
                AuxFlip = true;
            end
            Aux = [1 3];
        elseif Ratio(1) >= Ratio(2) && Ratio(4) >= Ratio(3)
            % Signals of pol. V and H may be switched:
            if Ratio(4) >= Ratio(1)
                AuxFlip = true;
            end     
            Aux = [2 4];
        else
            Aux = [2 4];
        end
        % Updating control variables:
        Ratio(Aux) = [] ; Pos(Aux) = [];
        if AuxFlip
            % Signals of pol. V and H are switched:
            Rx = fliplr(Rx) ; Pos = flipud(Pos) ; Ratio = flipud(Ratio);
        end
    end
    
    % Defining initial position according to the corr. with largest ratio:
    [~,IPos] = max(Ratio); PosIn = Pos(IPos); PosFin = PosIn+size(Rx,1)-1;
    if PosFin > size(Tx,1)
        % Updating final position:
        PosFin = size(Tx,1);
        
        % Adjusting sequence Rx so as Tx and Rx are the same length:
        PosInRx = PosIn - size(Tx,1) + size(Rx,1);  Rx = Rx(PosInRx:end,:);
    end
    
    % Synchronized symbol (Tx) and bit (TxBits) sequences:
    Tx = Tx(PosIn:PosFin,:) ; TxBits = TxBits((PosIn-1)*m+1:PosFin*m,:);
end