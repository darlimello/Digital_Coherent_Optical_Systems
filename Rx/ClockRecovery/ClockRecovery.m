function [Out,varargout] = ClockRecovery(In,PSType,NSymb,ParamCR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Out,varargout] = ClockRecovery(In,PSType,NSymb,ParamCR)                %
%                                                                         %
%  This function performs clock recovery in signal 'In' using a DPLL      %
% structure consisting on an interpolator, a TED, a loop filter, and a NCO%
% Signal 'In' must be obtained at 2 Sa/s. The function uses the Gardner   %
% TED for NRZ pulses, and the MG TED for Nyquist shaped pulses. The length%
% of the output sequence is limited to the number of transmitted symbols  %
% 'NSymb' times the number of samples per symbol at the input signal 'In' %
% (2 Sa/s).                                                               %
% *Notes:-The clock recovery is not done in the first two samples of 'In'.%
%        -This function is designed to be used before dynamic equalization%
%         However, in the presence of PMD, extra processing for tracking  %
%         polarization rotations is required;                             %
%                                                                         %
% Input:                                                                  %
%   In      = Signal obtained at 2 Sa/S of one pol. orientation in which  %
%             clock recovery will be performed. 'In' must be a column     %
%             vector;                                                     %
%   PSType  = Type of pulse shaping filter: 'NRZ' or 'Nyquist';           %
%   ParamCR = Struct that specifies parameters of the clock recovery:     %
%      - ParamCR.ki: Constant of the integral part of the loop filter;    %
%      - ParamCR.kp: Constant of the proportional part of the loop filter;%
%      - ParamCR.DPLLVarOut: Flag to enable (true) or disable (false)     %
%            internal variables of the DPLL (i.e., the output of the loop %
%            filter 'Wk', the NCO output 'Etamn', the fractional interval %
%            'mun', the base point 'mn', and TED output 'ek') as output of%
%            the function;                                                %
%                                                                         %
% Output:                                                                 %
%   xOut      = Signal obtained after clock recovery (column vector);     %
%   varargout = For 'ParamCR.DPLLVarOut = true', internal variables of the%
%               DPLL are also outputs of the function: 'Wk', 'Etamn',     %
%               'mun','mn', and 'ek';                                     % 
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Integral and proportional constants:
    ki = ParamCR.ki ; kp = ParamCR.kp;
        
    % Initializing variables:
    Etamn = 0.5 ; Wk = 1 ; LF_I = Wk ; mun = 0 ; n  = 3 ; mn = n;     
    Out = zeros(1,length(In)) ; Out(1:3) = In(1:3) ; LIn = length(In);
    if ParamCR.DPLLVarOut
        WkVec = [] ; EtamnVec = [] ; munVec = []; mnVec = [] ; ekVec = [] ;
    end
        
    % Clock recovery:
    while mn <= LIn       
        if mn == LIn;
            In(mn+1) = 0;
        end  
        
        % Cubic interpolator with Farrow architecture:
        Out(n) = In(mn-2)*(-1/6*mun^3 +   0*mun^2 + 1/6*mun + 0)+ ...
                 In(mn-1)*( 1/2*mun^3 + 1/2*mun^2 -   1*mun + 0)+ ...
                 In(mn)  *(-1/2*mun^3 -   1*mun^2 + 1/2*mun + 1)+ ...
                 In(mn+1)*( 1/6*mun^3 + 1/2*mun^2 + 1/3*mun + 0); 
                            
        % Generating the Ts-spaced timing error indication signal 'e_n':
        if mod(n,2)==1            
            switch PSType
                case 'Nyquist'
                    % TED - Nyquist pulses:
                    ek = abs(Out(n-1)).^2.*...
                         (abs(Out(n-2)).^2 - abs(Out(n)).^2);
                case 'NRZ'
                    % TED - NRZ pulses:
                    ek = real(conj(Out(n-1)).*(Out(n) - Out(n-2)));
            end         
            
            % Loop Filter:
            LF_I = ki*ek + LF_I; LF_P = kp*ek; Wk = LF_P + LF_I; 
            
            % Parameters of the DPLL:
            if ParamCR.DPLLVarOut
               WkVec   = [WkVec Wk]   ; EtamnVec = [EtamnVec Etamn];
               munVec  = [munVec mun] ; mnVec    = [mnVec mn];
               ekVec   = [ekVec ek];
            end         
        end
        
        % NCO - Base point 'mk' and fractional interval 'mu_n':
        if -1 < (Etamn - Wk) && (Etamn - Wk) < 0
            mn = mn + 1;            
        elseif (Etamn - Wk) >= 0
            mn = mn + 2;
        end  
        Etamn = mod(Etamn - Wk,1) ; mun  = Etamn/Wk; 
                
        % Updating the temporal index 'n':  
        n = n + 1;            
    end      

    % Limiting the length of the output to NSymb*2:
    if NSymb*2 < length(Out)
        Out = Out(1:NSymb*2).';
    else
        Out = Out.';
    end
    
    % Parameters of the DPLL:
    if ParamCR.DPLLVarOut
        varargout{1} = WkVec.'  ; varargout{2} = EtamnVec.';
        varargout{3} = munVec.' ; varargout{4} = mnVec.'   ;
        varargout{5} = ekVec.'  ;
    end    
end