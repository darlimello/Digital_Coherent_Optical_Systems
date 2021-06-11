function [AIR] = AIR_SDAWGN(SNRdB,ModFormat,AIREval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIRAWGN [AIR] = AIR_SDAWGN(SNRdB,ModFormat,AIREval)                     %
%                                                                         %
%   This function calculates approximations of the theoretical AIRs for   %
% schemes with SD-SW and SD-BW decoders under the circularly symmetric    %
% Gaussian channel. The approximation is given by a 10-point Gauss-Hermite%
% Quadrature. The supported modulation formats are 'QPSK' and '16QAM'. For%
% SD-SW decoders, the AIR is based on the MI, while for SD-BW decoders,   %
% the AIR is based on the GMI. For SD-BW decoders, Gray labeling is       %
% assumed. The binary labels are equal to the ones used in the function   %
% 'SymbolGeneration'.                                                     %
%                                                                         %
% Input:                                                                  %
%   SNRdB     = SNR in dB;                                                %
%   ModFormat = Modulation format: 'QPSK' or '16QAM';                     %
%   AIREval   = Defines the AIR to the evaluated: 'SD-SW' or 'SD-BW';     % 
%                                                                         %
% Output:                                                                 %
%   AIR = Approximation of the theoretical AIR (in bit/symbol) for SD-SW  %
%         or SD-BW decoders (depending on the variable 'AIREval') under a %
%         circularly symmetric Gaussian channel.                          %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Defining parameters according to the modulation format:
    switch ModFormat  
        case 'QPSK'
            % Modulation order and number of bits per symbol:
            M = 4 ; m = log2(M);
            % Constellation symbols:
            x = [-1+1i; 1+1i; -1-1i; 1-1i] ; x = x/sqrt(2);            
            % Bit-mapping:
            bmap = [0 1; 0 0; 1 1; 1 0];
        case '16QAM'
            % Modulation order and number of bits per symbol:
            M = 16 ; m = log2(M);
            % Constellation symbols:          
            x = [3+3i; 1+3i; -1+3i; -3+3i; 3+1i; 1+1i; -1+1i; -3+1i;...
                 3-1i; 1-1i; -1-1i; -3-1i; 3-3i; 1-3i; -1-3i; -3-3i];
            x = x/sqrt(10);            
            % Bit-mapping:
            bmap = [0 0 0 0; 0 1 0 0; 0 1 0 1; 0 0 0 1; 1 0 0 0;...
                    1 1 0 0; 1 1 0 1; 1 0 0 1; 1 0 1 0; 1 1 1 0;...
                    1 1 1 1; 1 0 1 1; 0 0 1 0; 0 1 1 0; 0 1 1 1; 0 0 1 1];
    end

    % Constants of the 10-point Gauss Hermite Quadrature:
    Zeta  = [-3.43615911883773760; -2.53273167423278980; ...
             -1.75668364929988177; -1.03661082978951365; ...
             -0.34290132722370461;  0.34290132722370461; ...
              1.03661082978951365;  1.75668364929988177; ...
              2.53273167423278980;  3.43615911883773760];	
    Gamma = [0.76404328552326206e-5; 0.13436457467812327e-2; ...
             0.33874394455481063e-1; 0.24013861108231469;    ...
             0.61086263373532580;    0.61086263373532580;    ...
             0.24013861108231469;    0.33874394455481063e-1; ...
             0.13436457467812327e-2; 0.76404328552326206e-5];

     % Calculating the standard deviation:
     SNR = 10.^(SNRdB/10) ; Stdev = sqrt(mean(abs(x).^2)./(2*SNR));

     switch AIREval
     case 'SD-SW'
        % MI Calculation (Gauss-Hermite quadrature):
        AIR = zeros(numel(SNRdB),1);
        for n = 1:numel(Stdev)
        Sum = 0; 
        for i = 1:M
            for l_1 = 1:numel(Zeta)
                for l_2 = 1:numel(Zeta)
                    num = sum(exp(-(abs(x(i)-x).^2 + 2*sqrt(2)*Stdev(n)*...
                        real((Zeta(l_1)+1i*Zeta(l_2)).*(x(i)-x)))/...
                        (2*Stdev(n)^2)));                
                    Sum = Sum + Gamma(l_1)*Gamma(l_2)*log2(num);
                end
            end
        end
        AIR(n) = m - 1/(M*pi)*Sum;
        end 
     case 'SD-BW'
        % GMI Calculation (Gauss-Hermite quadrature):
        AIR = zeros(numel(SNRdB),1);
        for n = 1:numel(Stdev)
        Sum = 0; 
        for k = 1:m
             for b = 0:1
                xSet = x(bmap(:,m+1-k) == b);
                for i = 1:M/2
                    for l_1 = 1:numel(Zeta) 
                        for l_2 = 1:numel(Zeta)
                            % Numerator
                            num = sum(exp(-(abs(xSet(i)-x).^2 +...
                                2*sqrt(2)*Stdev(n)*real((Zeta(l_1)+1i*...
                                Zeta(l_2)).*(xSet(i)-x)))/(2*Stdev(n)^2)));
                            % Denominator
                            den = sum(exp(-(abs(xSet(i)-xSet).^2 +...
                                2*sqrt(2)*Stdev(n)*real((Zeta(l_1)+1i*...
                                Zeta(l_2)).*(xSet(i)-xSet)))/...
                                (2*Stdev(n)^2)));
                            Sum = Sum + Gamma(l_1)*Gamma(l_2)*...
                                log2(num/den);               
                        end
                    end
                end
             end
         end
         AIR(n) = m - 1/(M*pi)*Sum;
        end
     end  
end