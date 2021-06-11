function [AOutput] = CDInsertion(AInput,SpS,Rs,D,CLambda,L,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CDInsertion [AOutput] = CDInsertion(AInput,SpS,Rs,D,CLambda,L,NPol)   %
%                                                                         %
%   This function simulates the insertion of chromatic dispersion (CD) in %
% the signal 'AInput'.                                                    %
%                                                                         %
% Input:                                                                  %
%  AInput = Input Signal. For transmission in single pol. orientation     %
%           'AInput' must be a column vector. For transmission with pol.  %
%           multiplexing, 'AInput' must be a matrix with two columns,     %
%           where each column corresponds to the signal of a pol.         %
%           orientation (V and H pol. orientations);                      %
%  SpS     = Number of samples per symbol in the input signal 'AInput';   %
%  Rs      = Symbol rate in [symbols/s];                                  %
%  D       = Dispersion parameter in [ps/(nm x km)]                       %
%  CLambda = Central lambda in [m]                                        %
%  L       = Fiber length in [m]                                          %
%  NPol    = Number of polarizations used;                                %
%                                                                         %
% Output:                                                                 %
%   AOutput = Output signal after CD insertion. 'AOutput' is arranged in  %
%             columns in the same way as 'AInput';                        %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Constants:
    c = 299792458;

    % Dispersion:
    D = D*1e-12/(1e-9*1e3);
    
    % Frequency vector:
    w = 2*pi*(-1/2:1/size(AInput,1):1/2-1/size(AInput,1)).'*SpS*Rs;
    
    % Calculating the CD frequency response:
    G = exp(1i*((D*CLambda^2)/(4*pi*c))*L*w.^2);

    % Inserting CD to the transmitted signal:
    AOutput(:,1) = ifft(ifftshift(G.*fftshift(fft(AInput(:,1)))));
    
    % In the case of pol-mux:
    if NPol == 2
    % Inserting CD to the transmitted signal:
        AOutput(:,2) = ifft(ifftshift(G.*fftshift(fft(AInput(:,2)))));
    end
end  