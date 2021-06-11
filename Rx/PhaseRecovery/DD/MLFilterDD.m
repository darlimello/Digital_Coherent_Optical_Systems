function [wML] = MLFilterDD(Delta_nu,Rs,OSNRdB,Es,NPol,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLFilterDD [wML] = MLFilterDD(Delta_nu,Rs,OSNRdB,Es,NPol,N)             %
%                                                                         % 
%  This function calculates the maximum likelihood (ML) filter for the    %
% Decision-Directed algorithm, which is calculated as                     %
%                         wML = (1^T C^-1)^T                              %
% where (.)^T indicates transpose and C is a covariance matrix, which     %
% depends on the signal-to-noise ratio and on the phase noise magnitude.  %
%                                                                         %
% Input:                                                                  %
%   Delta_nu = Sum of the transmitter and local oscillator laser          %
%              linewidths in Hz;                                          %
%   Rs       = Symbol rate in Symbols/second;                             %
%   OSNRdB   = Channel OSNR in dB;                                        %
%   Es       = Symbol energy (per po. orientation) in W;                  %
%   NPol     = Number of pol. orientations used;                          %
%   N        = Number of symbols used in the Decision-Directed algorithm  %
%              for phase noise estimation;                                %
%                                                                         %
% Output:                                                                 %
%   wML = Maximum likelihood filter to be used in the Decision-Directed   %
%         algorithm for phase noise estimation. wML is a column vector;   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Symbol period and phase noise variance:
    Ts = 1/Rs ; Sigma_DeltaTheta = 2*pi*Delta_nu*Ts;

    % Additive noise variance:
    SNRLin=10^(OSNRdB/10)*(2*12.5e9)/(NPol*Rs) ; Sigma_eta = Es/(2*SNRLin);

    % K matrix:
    K = zeros(N);
    for i = 0:N-1
        for j = 0:N-1
            K(i+1,j+1) = min(i,j);
        end
    end

    % Identity matrix:
    I = eye(N);

    % Obtaining the covariance matrix:
    C = Es^2*Sigma_DeltaTheta*K + Es*Sigma_eta*I;

    % Filter coefficients:
    wML = (ones(N,1)'/(C)).' ; wML = wML/max(wML);
end