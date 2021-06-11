function [E] = Laser(ParamLaser,SpS,Rs,NSymb,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASER [E] = Laser(ParamLaser,SpS,Rs,NSymb,NPol)                         %
%                                                                         %
%  This function simulates a laser as a continuous wave optical source,   % 
% according to the parameters defined in 'ParamLaser'. The electric field %
% 'E' of the optical signal has total power defined by 'ParamLaser.Pcw',  %
% for both the scenarios with one or two polarization orientations. Phase %
% noise is inserted into the electric field 'E' if the laser linewidth    %
% defined by 'ParamLaser.Linewidth' is different from 0 Hz.               %
%                                                                         %
% Input:                                                                  %
%   SpS        = Number of samples per symbol at the oversampled signal   %
%                (e.g., at the signal to be transmitted);                 %
%   Rs         = Symbol rate in symbols/second;                           %
%   NSymb      = Number of transmitted symbols;                           %
%   NPol       = Number of polarization orientations used;                %
%   ParamLaser = Struct that specifies parameters of the laser:           %
%     -ParamLaser.Pcw: Total power in dBm;                                %
%     -ParamLaser.Linewidth: Laser linewidth in Hz. The default value is  %
%              0 Hz (no phase noise). If the laser linewidth is not 0 Hz, %
%              phase noise is inserted into the optical carrier*;         %
%     *Note: 'ParamLaser.Linewidth' should be defined only if the required%
%            value differs from 0 Hz;                                     %
%                                                                         %
% Output:                                                                 %
%  E = Electric field of the optical signal. 'E' is a column vector if    %
%     NPol = 1, or a matrix with two column-oriented vectors if NPol = 2, %
%     where each column has the signal of a pol. orientation (V and H pol.%
%     orientations);                                                      %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Power of the continuous wave (Pcw) in dBm:
    Pcw = ParamLaser.Pcw;
    
    % Calculating the linear power of the continuous wave:
    PcwLinear = 1e-3*10^(Pcw/10);
    
    % Generating the electric filed of the optical signal:
    if NPol == 1      
        E = ones(SpS*NSymb,1)*sqrt(PcwLinear);         
    elseif NPol == 2        
        E = ones(SpS*NSymb,2)*sqrt(PcwLinear/2);         
    else        
        error('The possible number of polarizations used must be 1 or 2');        
    end
    
    % Laser Linewidth. Note: If the laser linewidth is not 0 Hz, phase 
    % noise is inserted into the optical carrier:
    if isfield(ParamLaser,'Linewidth')
        if ParamLaser.Linewidth ~= 0
            % Period between samples at the (oversampled) transmit. signal:
            T = 1/(SpS*Rs);

            % Calculating the phase noise:
            Var         = 2*pi*ParamLaser.Linewidth*T ; 
            Delta_theta = sqrt(Var)*randn(size(E,1),1);
            Theta       = cumsum(Delta_theta);

            % Adding phase noise to the optical signal:
            E = E.*repmat(exp(1i*Theta),1,size(E,2));
        end   
    end
end