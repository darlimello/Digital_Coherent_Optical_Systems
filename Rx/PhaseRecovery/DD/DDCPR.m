function [s_Decided,varargout] = DDCPR(z,ModFormat,Delta_nu,Rs,OSNRdB,Es...
    ,NPol,ParamDD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDCPR [s_Decided,varargout] = DDCPR(z,ModFormat,Delta_nu,Rs,OSNRdB,Es...%
%                                    ,NPol,ParamDD)                       %
%                                                                         %
%  This function performs phase recovery in signal 'z' with modulation    %
% format defined by 'ModFormat' using the Decision-Directed algorithm. The%
% algorithm uses 'ParamDD.N' previous symbol decisions to estimate the    %
% phase noise. In addition, a maximum likelihood (ML) filter is used prior%
% to generate phase estimates.                                            %
%                                                                         %
% Input:                                                                  %
%   z         = Signal in which phase recovery will be performed. For     %
%               transmission in single pol. orientation, 'z' must be a    %
%               column vector. For transmission with pol. multiplexing,'z'%
%               must be a matrix with two column vectors, where each      %
%               column vector corresponds to the signal of one pol.       %
%               orientation. Signal 'z' must be obtained at 1 sample per  %
%             symbol. Signal 'z' must also be normalized to unitary power;%
%   ModFormat = Modulation format of the signal 'z';                      %
%   Delta_nu  = Sum of the transmitter and local oscillator laser         %
%               linewidths in Hz;                                         %
%   Rs        = Symbol rate in symbols/second;                            %
%   OSNRdB    = Channel OSNR in dB;                                       %
%   Es        = Symbol energy (per pol. orientation) in W;                %
%   NPol      = Number of pol. orientations used;                         %
%   ParamDD   = Struct that specifies parameters of the decision-directed %
%               phase recovery algorithm:                                 %
%       -ParamDD.N: Number of past symbols used in the algorithm for phase%
%        noise estimation;                                                %
%       -ParamDD.OutSymb: Flag to enable ('true') or disable ('false') the%
%        symbols prior to decision as output;                             %
%       -ParamDD.PEstimate: Flag to enable ('true') or disable ('false')  %
%        the estimated phase noise as an output of the function;          %
%                                                                         %
% Output:                                                                 %
%   s_Decided = Sequence of (decided) symbols produced after compensating %
%               for the phase noise present on 'z';                       %
%   varargout: (1) When the flag 'ParamDD.OutSymb' is 'true', the symbols %
%         obtained prior to decision, 'v', are an output of the function; %
%              (2) When the flag 'ParamDD.PEstimate' is 'true', the       %
%         estimated phase noise is also an output of the function;        %
% Note:-Both the flags 'ParamDD.OutSymb' and 'ParamDD.PEstimate' can be   %
%    set 'true' at the same execution. In this case the function has 3    %
%    outputs:                                                             %
%                   [s_Decided,v,PEstimate] = DDCPR(....)                 %
%    - If only the flag 'ParamDD.OutSymb' is 'true',the outputs are:      %
%                      [s_Decided,v] = DDCPR(....)                        %
%    - If only the flag 'PEstimate' is 'true', the outputs are:           %
%                   [s_Decided,PEstimate] = DDCPR(....)                   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Vectors of phase estimates and output symbols:
    Theta = zeros(size(z,1),NPol) ; s_Decided = zeros(size(z,1),NPol);
        
    % Initializing the buffer of symbols used for phase recovery:
    N = ParamDD.N ; b = ones(N,NPol); 
    
    % ML filter:
    wML = MLFilterDD(Delta_nu,Rs,OSNRdB,Es,NPol,N);    
    
    % Phase noise estimates:
    for i = 1:length(z)               
        % Obtaining the phase of the previous received symbol:
        Theta(i,:) = angle(wML'*b);
        
        % Compensating for the phase of the received symbols:
        v = z(i,:).*exp(-1i*Theta(i,:));

        % Deciding the symbols after phase recovery:
        s_Decided(i,:) = Decision(v,ModFormat,false);

        % Updating the buffer b = z*s_Decided:
        b = circshift(b,1,1);
        b(1,:) = z(i,:).*conj(s_Decided(i,:))...
            ./abs(z(i,:).*conj(s_Decided(i,:)));
    end   
    
    % Symbols prior to decision:
    AuxOut = 1;
    if ParamDD.OutSymb
        varargout{AuxOut} = z.*exp(-1i*Theta) ; AuxOut = AuxOut + 1;
    end
    
    % Estimated phase noise as an output of the function:
    if ParamDD.PEstimate
       varargout{AuxOut} = unwrap(Theta); 
    end
end