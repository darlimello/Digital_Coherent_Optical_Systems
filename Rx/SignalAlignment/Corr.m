function [Pos,Ratio] = Corr(Tx,Rx,AbsCorr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corr [Pos,Ratio] = Corr(Tx,Rx,AbsCorr)                                  %
%                                                                         %
%   This function calculates cross-correlation between the transmitted    %
% sequence, 'Tx', and the received sequence, 'Rx', and produces 2 outputs.%
% The first one is the position 'Pos' of the maximum value of             %
% cross-correlation between 'Tx' and 'Rx'. The second one is the ratio    %
% between the maximum and the mean value of the cross-correlation, 'Ratio'%
%                                                                         %
% Input:                                                                  %
%  Tx      = Sequence of transmitted symbols normalized to unitary power*.%
%  Rx      = Sequence of received symbols normalized to unitary power*.   %
%  AbsCorr = Flag that determines if the cross-correlation is calculated  %
%            using only the magnitude ('true') or the complex values      %
%            ('false') of the symbols of sequences 'Tx' and 'Rx'. For     %
%            QPSK signals, 'AbsCorr' must be always 'false'. For 16-QAM   %
%            signals, a better synchronization is achieved when 'AbsCorr' %
%            is 'true';                                                   %
%  *Note: Sequences 'Tx' and 'Rx' are column vectors, corresponding to the%
%         signals of one pol. orientation;                                %
%                                                                         %
% Output:                                                                 %
%   Pos   = Position of the maximum value of the cross-correlation between%
%           'Tx' and 'Rx'.                                                %
%   Ratio = Ratio between the maximum and the mean value of the cross-    %
%           correlation between 'Tx' and 'Rx';                            %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Cross-correlation calculation:
    if ~AbsCorr
        Corr = xcorr(Tx-mean(Tx),Rx-mean(Rx));
    else
        Corr = xcorr(abs(Tx)-mean(abs(Tx)),abs(Rx)-mean(abs(Rx)));
    end
    
    % Finding the position of the maximum of the cross-correlation:
    Corr = Corr(size(Tx,1):end) ; [~,Pos] = max(abs(Corr));
    
    % Ratio between maximum and mean values:
    Ratio = 10*log10(max(abs(Corr))/mean(abs(Corr)));
end