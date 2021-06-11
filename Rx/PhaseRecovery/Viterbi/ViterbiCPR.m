function [v,varargout] = ViterbiCPR(z,Delta_nu,Rs,OSNRdB,Es,NPol,M,...
    ParamViterbi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VITERBICPR [v,varargout] = ViterbiCPR(z,Delta_nu,Rs,OSNRdB,Es,NPol,M,...%
%                                 ParamViterbi)                           %
%                                                                         %
%  This function performs phase recovery of signals with M-PSK modulation %
% formats using the Viterbi & Viterbi algorithm. The signal 'z' is first  %
% raised to M-th power and then ML filtering is performed before phase    %
% estimate. The phase estimates are applied to a phase unwrapper and,     %
% finally, the phase noise is compensated.                                %
%                                                                         %
% Input:                                                                  %
%   z            = Signal in which phase recovery will be performed. For  %
%                  transmission in single pol. orientation, 'z' must be a %
%                  column vector. For transmittion with pol. multiplexing,%
%                  'z' must be a matrix with two column-oriented vectors, %
%                  where each column vector corresponds to the signal of  %
%                  one pol. orientation. Signal 'z' must be obtained at 1 %
%                  sample per symbol. Signal 'z' must also be normalized  %
%                  to unitary power;                                      %  
%   Delta_nu     = Sum of the transmitter and local oscillator laser      %
%                  linewidths in Hz;                                      %
%   Rs           = Symbol rate in symbols/second;                         %
%   OSNRdB       = Channel OSNR in dB;                                    %
%   Es           = Symbol energy (per pol. orientation) in W;             %
%   NPol         = Number of pol. orientations used;                      %
%   M            = Modulation order of the M-PSK modulation format;       %
%   ParamViterbi = Struct that specifies parameters of the Viterbi        %
%                  & Viterbi algorithm:                                   %
%       -ParamViterbi.N: Number of past and future symbols used in the    %
%        Viterbi & Viterbi algorithm for phase noise estimate. The block  %
%        length is then 'L = 2*N+1';                                      %
%       -ParamViterbi.PEstimate: Flag to enable ('true') or disable       %
%        ('false') the estimated phase noise as an output of the function;%
%                                                                         %
% Output:                                                                 %
%   v = Signal produced after compensating for the phase noise present on %
%       'z';                                                              %
%   varargout: When the flag 'ParamViterbi.PEstimate' is 'true', the      %
%              estimated phase noise is also an output of the function;   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Block length for the phase estimation:
    N = ParamViterbi.N ; L = 2*N+1;
    
    % ML filter:
    wML = MLFilterViterbi(M,Delta_nu,Rs,OSNRdB,Es,NPol,N);

    % Initializing the vector of phase estimates:
    ThetaML = zeros(size(z,1),NPol);
    
    % Phase noise estimation:
    for Pol = 1:NPol
        % Input blocks:
        zBlocks = [zeros(floor(L/2),1) ; z(:,Pol) ; zeros(floor(L/2),1)];
        zBlocks=convmtx(zBlocks.',L); zBlocks=flipud(zBlocks(:,L:end-L+1));  

        % Note that each column of zBlocks is used independently to
        % generate phase estimates:
        ThetaML(:,Pol) = (1/M)*angle(wML.'*(zBlocks.^M))-pi/M;
    end
    clearvars zBlocks;
    
    % Vector of phase estimates after phase unwrapping:
    ThetaPU = zeros(size(ThetaML,1),NPol);
    
    % Initial 'previous phase' for unwrapping operation:
    ThetaPrev = zeros(1,NPol);
    
    % Phase unwrapping:
    for i = 1:size(ThetaML,1)
        % Phase unwrapper:
        n            = floor(1/2 + (ThetaPrev - ThetaML(i,:))/(2*pi/M));
        ThetaPU(i,:) = ThetaML(i,:) + n*(2*pi/M); ThetaPrev = ThetaPU(i,:);       
    end  
    
    % Phase noise compensation:
    v = z.*exp(-1i*ThetaPU);     
    
    % Estimated phase noise as an output of the function:
    if ParamViterbi.PEstimate
       varargout{1} = ThetaPU; 
    end
end