function [wML] = MLFilterViterbi(M,Delta_nu,Rs,OSNRdB,Es,NPol,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLFILTER [wML] = MLFilterViterbi(M,Delta_nu,Rs,OSNRdB,Es,NPol,N)        %
%                                                                         %
%  This function calculates the maximum likelihood (ML) filter for the    %
% Viterbi & Viterbi algorithm, which is calculated as                     %
%                           wML = (1^T C^-1)^T                            %
% where (.)^T indicates transpose and C is a covariance matrix, which     %
% depends on the signal-to-noise ratio and on the phase noise magnitude.  %
%                                                                         %
% Input:                                                                  %
%   M        = Modulation order of the M-PSK modulation format;           %
%   Delta_nu = Sum of transm. and local oscillator laser linewidths in Hz;%
%   Rs       = Symbol rate in symbols/second;                             %
%   OSNRdB   = Channel OSNR in dB;                                        %
%   Es       = Symbol energy (per pol. orientation) in W;                 %
%   NPol     = Number of pol. orientations used;                          %
%   N        = Number of past and future symbols used in the Viterbi &    %
%              Viterbi algorithm for phase noise estimates. The block     %
%              length is then L = 2*N+1;                                  %
%                                                                         %
% Output:                                                                 %
%   wML = Maximum likelihood filter to be used in the Viterbi & Viterbi   %
%         algorithm for phase noise estimation. wML is a column vector;   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Block length for the Viterbi & Viterbi algorithm, symbol period, and 
    % phase noise variance:
    L = 2*N+1 ; Ts = 1/Rs ; Sigma_DeltaTheta2 = 2*pi*Delta_nu*Ts;

    % Additive noise variance:
    SNRLin=10^(OSNRdB/10)*(2*12.5e9)/(NPol*Rs) ; Sigma_eta2= Es/(2*SNRLin);

    % K matrix:
    KAux = zeros(N) ; K    = zeros(L);
    for i = 0:N
        for ii = 0:N
            KAux(i+1,ii+1) = min(i,ii);
        end
    end
    K(1:N+1,1:N+1) = rot90(KAux(1:N+1,1:N+1),2);
    K(N+1:L,N+1:L) = KAux(1:N+1,1:N+1);

    % Identity matrix:
    I = eye(L);

    % Obtaining the covariance matrix:
    C = Es^M*M^2*Sigma_DeltaTheta2*K + Es^(M-1)*M^2*Sigma_eta2*I;

    % Filter coefficients:
    wML = (ones(L,1)'/(C)).' ; wML = wML/max(wML);
end