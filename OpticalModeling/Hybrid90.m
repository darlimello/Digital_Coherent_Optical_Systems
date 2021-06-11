function [E1,E2,E3,E4] = Hybrid90(Er,ELo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYBRID90 [E1,E2,E3,E4] = Hybrid90(Er,ELo)                               %
%                                                                         %
%   This function simulates a 90 degree hybrid (with ideal components).   %
%                                                                         %
% Input:                                                                  %
%   Er  = Received optical signal in one pol. orientation (column vector);%
%   ELo = Local oscillator signal in one pol. orientation (column vector);%
%                                                                         %
% Output:                                                                 %
%   E1, E2, E3, and E4 = 90 degree hybrid output signals;                 %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3-dB coupler transfer function:
    Hc = (1/sqrt(2))*[1 1; 1 -1];
    
    % ECouplerTL - Signal at the output of the top-left 3-dB coupler at the
    % 90 degree hybrid;
    ECouplerTL = (Hc*[Er.' ; zeros(1,length(Er))]).';
    
    % ECouplerBL - Signal at the output of the bottom-left 3-dB coupler at 
    % the 90 degree hybrid;
    ECouplerBL = (Hc*[ELo.' ; zeros(1,length(Er))]).';
    
    % ECouplerTR - Signal at the output of the top-right 3-dB coupler at 
    % the 90 degree hybrid;
    ECouplerTR = (Hc*[ECouplerTL(:,1).' ; ECouplerBL(:,1).']).';
    
    % ECouplerBR - Signal at the output of the bottom-right 3-dB coupler at
    % the 90 degree hybrid;
    ECouplerBR = (Hc*[ECouplerTL(:,2).';ECouplerBL(:,2).'*exp(1i*pi/2)]).';
    
    % Output signals:
    E1 = ECouplerTR(:,1); E2 = ECouplerTR(:,2); E3 = ECouplerBR(:,1);
    E4 = ECouplerBR(:,2);
end