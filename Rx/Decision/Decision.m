function [Decided] = Decision(r,ModFormat,BitsOutput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECISION [Decided] = Decision(r,ModFormat,BitsOutput)                   %
%                                                                         %
%  This function performs decisions on each symbol of sequence 'r'        %
% according to the modulation format 'ModFormat'. The function generates  %
% as output the sequence of decided symbols or the sequence of bits       %
% corresponding to the decided symbols (assuming the same binary labeling %
% of function 'SymbolGeneration'), depending on the flag 'BitsOutput'.    %
% This function does not perform decisions in the case of differential    %
% encoding.                                                               %
%                                                                         %
% Input:                                                                  %
%   r          = Sequence of symbols received (one pol. orientation -     % 
%                column vector) normalized to unitary power;              %
%   ModFormat  = Modulation format: 'QPSK' or '16-QAM';                   %
%   BitsOutput = If 'BitsOutput = false', the function generates as output%
%                the sequence of decided symbols. If 'BitsOutput = true', %
%                the function generates as output the sequence of bits    %
%                corresponding to the decided symbols;                    %
%                                                                         %
% Output:                                                                 %
%   Decided = Sequence of decided symbols (if 'BitsOutput = false') or    %
%             the sequence of bits corresponding to the decided symbols   %
%             (if 'BitsOutput = true'). In both cases, 'Decided' is a     %
%             column vector.                                              %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Decision:
    switch ModFormat
        case 'QPSK'
            % Decision regions for the in-phase components:
            R1 = real(r) >= 0; R2 = real(r) < 0;
            % Applying the decision regions to the imaginary axis:
            R3 = imag(r) >= 0; R4 = imag(r) < 0;
        case '16QAM'
            % Applying the decision regions to the real axis:
            R1 = real(r) >= 2/sqrt(10); R2 = real(r) >= 0;
            R3 = real(r)  < 0         ; R4 = real(r) <=-2/sqrt(10);
            % Applying the decision regions to the imaginary axis:
            R5 = imag(r) >= 2/sqrt(10); R6 = imag(r) >= 0;
            R7 = imag(r) < 0          ; R8 = imag(r) <=-2/sqrt(10);
        otherwise
            error('The Supported modulation formats are QPSK and 16-QAM;');
    end
            
    if BitsOutput      
        % Binary labeling:
        switch ModFormat                    
            case 'QPSK'
                ModBits = 2;
                Decided = NaN(length(r),2);
                % Assigning the bits based on the mapping done in the 
                % transmitter:
                Decided(R1,2) = zeros(1,sum(R1));
                Decided(R2,2) = ones(1,sum(R2));
                Decided(R3,1) = zeros(1,sum(R3));
                Decided(R4,1) = ones(1,sum(R4));
            case '16QAM'
                ModBits = 4;
                Decided = NaN(length(r),4);
                % Assigning the bits based on the mapping done in the 
                % transmitter:
                Decided(R1,[2 4])     = repmat([0 0],sum(R1),1);
                Decided(R2&~R1,[2 4]) = repmat([1 0],sum(R2&~R1),1);
                Decided(R3&~R4,[2 4]) = repmat([1 1],sum(R3&~R4),1);
                Decided(R4,[2 4])     = repmat([0 1],sum(R4),1);
                Decided(R5,[1 3])     = repmat([0 0],sum(R5),1);
                Decided(R6&~R5,[1 3]) = repmat([1 0],sum(R6&~R5),1);
                Decided(R7&~R8,[1 3]) = repmat([1 1],sum(R7&~R8),1);
                Decided(R8,[1 3])     = repmat([0 1],sum(R8),1);           
        end 
        % Obtaining the decided bits as a column vector:
        Decided = reshape(Decided',1,length(r)*ModBits)';
    else
        % Decided Symbols:
        Decided = zeros(size(r));
        switch ModFormat            
            case 'QPSK'                
                % Assigning the bits based on the mapping done in the 
                % transmitter:
                Decided(R1) = Decided(R1) + (1/sqrt(2));
                Decided(R2) = Decided(R2) - (1/sqrt(2));
                Decided(R3) = Decided(R3) + (1i*1/sqrt(2));
                Decided(R4) = Decided(R4) - (1i*1/sqrt(2));                  
            case '16QAM'                
                % Assigning the bits based on the mapping done in the 
                % transmitter:
                Decided(R1)       = Decided(R1)     + (3/sqrt(10));
                Decided(R2 & ~R1) = Decided(R2&~R1) + (1/sqrt(10));
                Decided(R3 & ~R4) = Decided(R3&~R4) - (1/sqrt(10));
                Decided(R4)       = Decided(R4)     - (3/sqrt(10));
                Decided(R5)       = Decided(R5)     + (3i/sqrt(10));
                Decided(R6 & ~R5) = Decided(R6&~R5) + (1i/sqrt(10));
                Decided(R7 & ~R8) = Decided(R7&~R8) - (1i/sqrt(10));
                Decided(R8)       = Decided(R8)     - (3i/sqrt(10));                                
        end
    end
end