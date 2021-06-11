function [rOut] = Deskew(rIn,SpSRx,Rs,NPol,N,ParamSkew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESKEW [rOut] = Deskew(rIn,SpSRx,Rs,NPol,N,ParamSkew)                   %
%                                                                         %
%  This function performs deskew in the signals 'rIn' using a Lagrange    %
% interpolator of order 'N'. The interpolator is implemented by a FIR     %
% filter of length 'N+1'. The temporal misalignment is compensated taking %
% into account the lowest temporal delay. The temporal delays of each     %
% of the input signal (in-phase and quadrature components) are specified  %
% in 'ParamSkew'.                                                         % 
%                                                                         %
% Input:                                                                  %
%   rIn       = Input signal in which deskew will be performed. For       %
%               transmission in single pol. orientation, 'In' must be a   %
%               matrix with two columns. The 1st and 2nd columns must     %
%               have the in-phase and quadrature components of the signal,%
%               respectively. For transmission with pol. multiplexing,    %
%               'In' must be a matrix with four columns. The 1st and the  %
%               2nd columns must have the in-phase and quadrature         %
%               components of the signal in V pol. orientation, while the %
%               3rd and 4th columns, the in-phase and quadrature          %
%               components of the signal in H pol. orientation,           %
%               respectively;                                             %
%   SpSRx     = Samples/symbol at the receiver (after the ADC);           %
%   Rs        = Symbol rate in Symbols/second;                            %
%   NPol      = Number of polarization orientations;                      %
%   N         = Order of the Lagrangean interpolation polynomial;         %
%   ParamSkew = Struct that specifies the temporal delay for each         %
%               component of the input signal:                            %
%                      ParamSkew.TauIV and ParamSkew.TauQV;               %
%                      ParamSkew.TauIH and ParamSkew.TauQH (if NPol == 2) %
%                                                                         %
% Output:                                                                 %
%   rOut = Signal after deskew. 'rOut' is organized in columns in the same%
%          way as 'rIn';                                                  %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Skew to be compensated using the TADC as reference:
    TADC=1/(SpSRx*Rs) ; Skew = [ParamSkew.TauIV/TADC ParamSkew.TauQV/TADC];
    if NPol == 2
        Skew = [Skew ParamSkew.TauIH/TADC ParamSkew.TauQH/TADC];
    end
    
    % Using the min skew as reference:
    Skew = Skew - min(Skew);

    % Integer and fractional part of the skew:
    nTADC = floor(Skew) ; muTADC = -(Skew-nTADC);  
    
    % Obtaining the FIR filter and interpolating the signals: 
    NTaps = N+1; % Number of filter taps;
    for i = 1:size(rIn,2)                        
        L = zeros(NTaps,1) ; Aux = 1;       
        
        % Obtaining the Lagrangean interpolator:
        for n = (0:N) - floor(mean(0:N)) + nTADC(i)
            m      = (0:N) - floor(mean(0:N)) + nTADC(i) ; m(m == n) = [];
            L(Aux) = prod((muTADC(i) - m)./(n - m))      ; Aux = Aux + 1;
        end
        
        % Interpolating the received signal (sIn):
        sAux = flipud(convmtx([zeros(1,floor(NTaps/2)) rIn(:,i).'...
            zeros(1,floor(NTaps/2))],NTaps));
        sAux = sAux(:,NTaps:end-(NTaps)+1) ; rOut(:,i) = (L.'*sAux).'; 
    end
end