function [r] = NoiseInsertion_OSNR(x,OSNRdB,SpS,NPol,Rs,BRef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISEINSERTION_OSNR [r] = NoiseInsertion_OSNR(x,OSNRdB,SpS,NPol,Rs,Bref)%
%                                                                         %
%   This function inserts additive white Gaussian noise (AWGN) in the     %
% transmitted signal 'x' (in single- or dual-polarization orientations) so%
% that an OSNR (in dB) of 'OSNRdB' is achieved.                           %
%                                                                         %
% Input:                                                                  %
%   x      = Transmitted signal. For transmission in single pol.          %
%            orientation, 'x' must be a column vector. For transmission   %
%            with pol. multiplexing, 'x' must be a matrix with two column,%
%            where each column vector corresponds to the signal of a pol. %
%            orientation (V and H pol. orientations);                     %
%   OSNRdB = OSNR in dB;                                                  %
%   SpS    = Number of samples per symbol in the input signal 'x';        %
%   NPol   = Number of polarizations to be used.                          %
%   Rs     = Symbol rate in symbols/second;                               %
%   BRef   = Reference bandwidth (in Hz) for OSNR measurement;            %
%                                                                         %
% Output:                                                                 %
%   r = Signal after noise insertion. 'r' is arranged in columns in the   %
%       same way as 'x';                                                  %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % OSNR/SNR in linear scale:
    OSNRLin = 10^(OSNRdB/10);
    SNRLin  = ((2*BRef)/(NPol*Rs))*OSNRLin;
    
    % AWGN standard deviation:
    StdDev = sqrt(mean(abs(x(:,1)).^2)*SpS/(2*SNRLin));
    
    % AWGN generation:
    n = StdDev*randn(length(x(:,1)),1)+1i*StdDev*randn(length(x(:,1)),1);
    
    % Inserting noise to the signal:    
    r(:,1) = x(:,1) + n;   
       
    % In the case of pol-mux:
    if NPol == 2
      % AWGN standard deviation:
      StdDev = sqrt(mean(abs(x(:,2)).^2)*SpS/(2*SNRLin));
        
      % AWGN generation:
      n = StdDev*randn(length(x(:,2)),1)+1i*StdDev*randn(length(x(:,2)),1);
            
      % Inserting noise to the signal:    
      r(:,2) = x(:,2) + n;
    end
end