function [Bits,x] = SymbolGeneration(ModFormat,NSymb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMBOLGENERATION [Bits,x] = SymbolGeneration(ModFormat,NSymb)           %
%                                                                         %
%  This function generates a sequence of bits 'Bits' with uniform         %
% distribution and, considering Gray mapping, generates a sequence of     %
% symbols 'x' according to the modulation format 'ModFormat'. The length  %
% of 'x' is given by 'NSymb'. The symbols of 'x' follow the unitary power %
% constellation associated with 'ModFormat'. For transmission with pol.   %
% multiplexing, this function must be called twice.                       %
%                                                                         %
% Input:                                                                  %
%   ModFormat  = Modulation format: 'QPSK' or '16QAM';                    %
%   NSymb      = Number of symbols to be transmitted;                     %
%                                                                         %
% Output:                                                                 %
%   Bits = Sequence of bits (column vector);                              %
%   x    = Sequence of symbols (column vector);                           %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch ModFormat        
        case {'QPSK'}
            ModBits = 2;
        case '16QAM'
            ModBits = 4;
        otherwise
            error('The supported modulation formats are QPSK and 16-QAM;');
    end

    % Generating the sequence of bits:
    NBits = NSymb*ModBits; Bits = randi([0,1],NBits,1);
        
    switch ModFormat
        case 'QPSK'
            % In-Phase and quadrature bits:
            BitsI = Bits(2:ModBits:NBits); BitsQ = Bits(1:ModBits:NBits);

            % QPSK modulation:
            xI = 1-2*BitsI; xQ = 1-2*BitsQ;
            
            % Normalized symbols:
            x = (1/sqrt(2))*(xI + xQ*1i);
        case '16QAM'
            % In-Phase and quadrature bits:
            BitsI1 = Bits(4:ModBits:NBits); BitsQ1 = Bits(3:ModBits:NBits);
            BitsI2 = Bits(2:ModBits:NBits); BitsQ2 = Bits(1:ModBits:NBits);

            % 16QAM modulation:
            % In-Phase:
            xI = ((~BitsI2 & ~BitsI1)*(+3) + ( BitsI2 & ~BitsI1)*(+1) + ...
                  ( BitsI2 &  BitsI1)*(-1) + (~BitsI2 &  BitsI1)*(-3)); 

            % Quadrature:
            xQ = ((~BitsQ2 & ~BitsQ1)*(+3) + ( BitsQ2 & ~BitsQ1)*(+1) + ...
                  ( BitsQ2 &  BitsQ1)*(-1) + (~BitsQ2 &  BitsQ1)*(-3));

            % Normalized symbols:
            x = (1/sqrt(10))*(xI + xQ*1i);
    end        
end