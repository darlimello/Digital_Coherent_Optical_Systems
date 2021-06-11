function [z,varargout] = FRecovery(y,Rs,FEstimate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOCOMPENSATION [z] = FRecovery(y,Rs,FEstimate)                          %
%                                                                         %
%  This function performs frequency recovery in the signal 'y'. Using the %
% spectal analysis of the signal 'y' raised to 4-th power an estimate     %
% 'Delta_f' of the carrier frequency offset is generated. Then, the       %
% frequency recovery is performed as                                      %
%                      z = y.*exp(-1i*k*2*pi*Delta_f*Ts)                  %
% where 'Ts' is the symbol period and 'k' is the temporal index of the    %
% symbols;                                                                %
%                                                                         %
% Input:                                                                  %
%   y         = Signal in which the frequency recovery will be performed. %
%               Signal 'y' must be obtained at 1 sample per symbol. For   %
%               transmission in single pol. orientation, 'y' must be a    %
%               column vector. For transmission with pol. multiplexing,   %
%               'y' must be a matrix with two column-oriented vectors,    %
%               where each column vector corresponds to the signal of one %
%               pol. orientation;                                         %
%   Rs        = Symbol rate in symbols/second;                            %
%   FEstimate = Flag to enable ('true') or disable ('false') the estimated%
%             carrier frequency offset value as an output of the function;%
%                                                                         %
% Output:                                                                 %
%   z         = Signal produced after compensating for the frequency      %
%               offset present on 'y';                                    %
%   varargout = Estimated carrier frequency offsert in Hz;                %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Frequency vector considering 1 sample per symbol and symbol period:
    f = (-1/2+1/length(y):1/length(y):1/2)*Rs ; Ts = 1/Rs;

    % Obtaining the absolute value of the spectrum of the signal^4: 
    SignalSpectrum = fftshift(abs(fft(y(:,1).^4)));

    % Obtaining the frequency offset:
    Delta_f = (1/4)*f(SignalSpectrum == max(SignalSpectrum));

    % Compensating for the frequency offset:
    k = repmat((0:length(y)-1).',1,size(y,2));
    z = y.*exp(-1i*2*pi*Delta_f*Ts*k);
    
    % Estimated carrier frequency offset as an output of the function:
    if FEstimate
       varargout{1} = Delta_f; 
    end    
end