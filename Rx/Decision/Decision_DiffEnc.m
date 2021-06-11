function [Decided] = Decision_DiffEnc(r,ModFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision_DiffEnc [Decided] = Decision_DiffEnc(r,ModFormat)              %
%                                                                         %
%  This function performs differential decoding on the received sequence  %
% 'r' according to the modulation format 'ModFormat'. As output, the      %
% function generates the binary sequence associated with the sequence 'r'.%
%                                                                         %
% Input:                                                                  %
%   r         = Received signal (one pol. orientation - column vector)    %
%               normalized to unitary power;                              %
%   ModFormat = Modulation format: 'QPSK' or '16-QAM';                    %
%                                                                         %
% Output:                                                                 %
%   Decided = Binary sequence associated with sequence 'r' (column vector)%
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    switch ModFormat        
        case 'QPSK'
            ModBits = 2;
        case '16QAM'
            ModBits = 4;
        otherwise
            error('The supported modulation formats are QPSK and 16-QAM;');
    end
    
    % Decision:
    QPrev = 0;
    Decided  = NaN(length(r),ModBits);
    switch ModFormat             
        case 'QPSK'
            % Decision regions for the in-phase component:
            R1 = real(r) >= 0 ; R2 = real(r) < 0;
            % Decision regions for the quadrature component:
            R3 = imag(r) >= 0 ; R4 = imag(r) < 0;

            % Defining the quadrant:
            Q(R1&R3) = 0 ; Q(R2&R3) = 1 ; Q(R2&R4) = 2 ; Q(R1&R4) = 3;
            
            % Binary sequence:
            for i = 1:length(r)
                QRx = (Q(i) - QPrev);
                % Bits defining the quadrant (and consequently the symbol):
                switch QRx
                    case 0
                        Decided(i,1) = 0; Decided(i,2) = 0;
                    case {1,-3}
                        Decided(i,1) = 0; Decided(i,2) = 1;
                    case {3,-1}
                        Decided(i,1) = 1; Decided(i,2) = 0;
                    case {2,-2}
                        Decided(i,1) = 1; Decided(i,2) = 1;
                end                            
                QPrev = Q(i);                         
            end          
        case '16QAM'
            % Applying the decision regions to the real axis:
            R1 = real(r) >= 2/sqrt(10); R2 = real(r) >= 0;
            R3 = real(r) < 0          ; R4 = real(r) <= -2/sqrt(10);
            % Applying the decision regions to the imaginary axis:
            R5 = imag(r) >= 2/sqrt(10); R6 = imag(r) >= 0;
            R7 = imag(r) < 0          ; R8 = imag(r) <= -2/sqrt(10);

            % Defining the quadrant:
            Q(R1&R5 | R1&R6&~R5 | ~R1&R2&R5 | ~R1&R2&R6&~R5) = 0;         
            Q(~R4&R3&R5 | ~R4&R3&R6&~R5 | R4&R5 | R4&R6&~R5) = 1;
            Q(R4&R7&~R8 | R4&R8 | ~R4&R3&R7&~R8 | ~R4&R3&R8) = 2;
            Q(~R1&R2&R7&~R8 | ~R1&R2&R8 | R1&R7&~R8 | R1&R8) = 3;

            % Defining the symbol inside the quadrants:
            S(~R1&R2&~R5&R6|~R1&R2&~R8&R7|~R4&R3&~R5&R6|~R4&R3&~R8&R7) = 0;
            S(~R1&R2&R5 | R4&R6&~R5 | ~R4&R3&R8 |  R1&R7&~R8)          = 1;
            S(R1&R6&~R5 | R4&R7&~R8 | ~R4&R3&R5 | ~R1&R2&R8)           = 2; 
            S(R1&R5     | R1&R8     | R4&R5     | R4&R8)               = 3;

            % Received binary sequence:
            for i = 1:length(r)
                % Bits defining the quadrant:
                QRx = (Q(i) - QPrev);
                switch QRx
                    case 0
                        Decided(i,1) = 0; Decided(i,2) = 0;
                    case {1,-3}
                        Decided(i,1) = 0; Decided(i,2) = 1;
                    case {3,-1}
                        Decided(i,1) = 1; Decided(i,2) = 0; 
                    case {2,-2}
                        Decided(i,1) = 1; Decided(i,2) = 1;                       
                end                            
                QPrev = Q(i);
                % Bits defining the symbol inside the quadrant:
                switch S(i)
                    case 0
                        Decided(i,3) = 0; Decided(i,4) = 0;  
                    case 1
                        Decided(i,3) = 0; Decided(i,4) = 1;  
                    case 2
                        Decided(i,3) = 1; Decided(i,4) = 0; 
                    case 3
                        Decided(i,3) = 1; Decided(i,4) = 1; 
                end                             
            end  
    end             
    % Obtaining the vector of bits:
    Decided = reshape(Decided',1,length(r)*ModBits)';    
end