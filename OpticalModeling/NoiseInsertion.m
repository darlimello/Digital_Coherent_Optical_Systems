function [r] = NoiseInsertion(x,ModBits,SNRb_dB,SpS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISEINSERTION [r] = NoiseInsertion(x,ModBits,SNRb_dB,SpS)              %
%                                                                         %
%  This function inserts additive white Gaussian noise (AWGN) in the      %
% transmitted signal 'x' (single-polarization orientation), so that an    %
% SNR per bit (in dB) 'SNRb_dB' is achieved.                              %
%                                                                         %
% Input:                                                                  %
%   x       = Transmitted signal (single pol. orientation - column vector)%
%   ModBits = Number of bits per symbol of the modulation format used in  %
%             in the signal 'x';                                          %
%   SNRb_dB = SNR per bit in dB                                           %
%   SpS     = Number of samples per symbol at the input signal 'x';       %
%                                                                         %
% Output:                                                                 %
%   r = Signal after noise insertion (column vector);                     %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SNR per bit in linear scale:
    SNRb_Lin = 10^(SNRb_dB/10); 
    
    % AWGN standard deviation:
    StdDev = sqrt(mean(abs(x).^2)*SpS/(2*ModBits*SNRb_Lin));
    
    % AWGN generation:
    n = StdDev*randn(length(x),1) + 1i*StdDev*randn(length(x),1);
    
    % Inserting noise to the signal:    
    r = x + n;            
end