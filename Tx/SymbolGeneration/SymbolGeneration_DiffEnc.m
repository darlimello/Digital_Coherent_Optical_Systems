function [Bits,x] = SymbolGeneration_DiffEnc(ModFormat,NSymb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SymbolGeneration_DiffEnc [Txbits,x] = SymbolGeneration_DiffEnc(...      %
%                               ModFormat,NSymb)                          %
%                                                                         %
%  This function generates a sequence of bits 'Bits' with uniform         %
% distribution and, considering differential encoding, generates a        %
% sequence of symbols 'x' according to the modulation format 'ModFormat'. %
% The length of 'x' is given by 'NSymb'. The symbols of 'x' follow the    %
% unitary power constellation associated with 'ModFormat'. For            %
% transmission with pol. multiplexing, this function must be called twice.%
%                                                                         %
% Input:                                                                  %
%   ModFormat  = Modulation format: 'QPSK' or '16QAM';                    %
%   NSymb      = Number of symbols to be transmitted;                     %
%                                                                         %
% Output:                                                                 %
%   Bits     = Bits to be transmitted (column vector);                    %
%   x        = Sequence of symbols (column vector);                       %
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

    % Generating bits to be transmitted:
    NBits = NSymb*ModBits; Bits = randi([0,1],NBits,1);
        
    switch ModFormat  
        case 'QPSK'
            % Input bits:
            QBits1 = Bits(1:ModBits:NBits); QBits2 = Bits(2:ModBits:NBits);

            % Defining the quadrant:
            Quad = ((~QBits1&~QBits2)*1            + ...
                    (~QBits1& QBits2)*exp(1i*pi/2) + ...
                    ( QBits1&~QBits2)*exp(3i*pi/2) + ...
                    ( QBits1& QBits2)*exp(1i*pi));

           % Initial quadrant:
           QuadPrev = 1;

           % Modulation:
           x = NaN(1,length(Quad));
           for i = 1:length(x)
               x(i)     = QuadPrev*Quad(i)*(1+1i); 
               QuadPrev = QuadPrev*Quad(i);
           end
           x = (1/sqrt(2))*x.'; 

        case '16QAM'
            % Bits that define the quadrant:
            QBits1 = Bits(1:ModBits:NBits); QBits2 = Bits(2:ModBits:NBits);

            % Bits that define the symbols inside the quadrants:
            InQuadBits1 = Bits(3:ModBits:NBits);
            InQuadBits2 = Bits(4:ModBits:NBits);              

            % Defining the quadrant:
            Quad = ((~QBits1&~QBits2)*1            + ...
                    (~QBits1& QBits2)*exp(1i*pi/2) + ...
                    ( QBits1&~QBits2)*exp(3i*pi/2) + ...
                    ( QBits1& QBits2)*exp(1i*pi)); 

            % Defining the symbol inside the quadrants:
            InQuadI      = 2*InQuadBits1 + 1; InQuadQ = 2*InQuadBits2 + 1;
            InQuadSymbol = InQuadI + 1i*InQuadQ;

            % Initial quadrant:
            QuadPrev = 1;

            % Modulation:
            x = NaN(1,length(InQuadBits2));  
            for i = 1:length(x)
               x(i)     = QuadPrev*Quad(i)*InQuadSymbol(i);
               QuadPrev = QuadPrev*Quad(i);
            end
            x = ((1/sqrt(10))*x).';
   end      
end