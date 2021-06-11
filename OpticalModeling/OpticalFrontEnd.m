function [Out] = OpticalFrontEnd(Er,ELo,R,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COHERENTDETECTORFE [Out] = OpticalFrontEnd(Er,ELo,R,NPol)               %
%                                                                         %
%   This function simulates a single- or dual-polarization optical        %
% front-end, considering ideal 90 degree hybrids.                         %
%                                                                         %
% Input:                                                                  %
%   Er   = Signal to be detected*;                                        %
%   ELo  = Signal produced by the local oscillator laser;                 %
%   R    = Photodetector responsivity;                                    %
%   NPol = Number of polarizations;                                       %
%  *Note = For transmission in single pol. orientation 'Er' and 'ELo' must%
%          be column vectors. For transmission with pol. multiplexing,    %
%          'Er' and 'ELo' must be matrices with with two columns, where   %
%          each column corresponds to the signal of a pol. orientation,   %
%          (V and H pol. orientations);                                   %
%                                                                         %
% Output:                                                                 %
%   Out = Electric currents produced at the front-end. For 'NPol = 1',    %
%         'Out' consists on the electric currents iIV and iQV (column     %
%         vectors). For 'NPol = 2', 'Out' consists on the electric        %
%         currents iIV, iQV, iIH, and iQH (column vectors);               %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 90º Hybrid:
    [E1,E2,E3,E4] = Hybrid90(Er(:,1),ELo(:,1));

    % Photodetection:
    i1 = R*(E1.*conj(E1)) ; i2 = R*(E2.*conj(E2)); % In-phase;
    i3 = R*(E3.*conj(E3)) ; i4 = R*(E4.*conj(E4)); % Quadrature;
    
    % Balanced photodetection:
    iI = i1 - i2 ; iQ  = i3 - i4;

    % Output signal:
    Out(:,1) = iI ; Out(:,2) = iQ;
    
    % In the case of pol-mux:
    if NPol == 2 
        % 90º Hybrid:
        [E1,E2,E3,E4] = Hybrid90(Er(:,2),ELo(:,2));    

        % Photodetection:
        i1 = R*(E1.*conj(E1)) ; i2 = R*(E2.*conj(E2)); % In-phase;
        i3 = R*(E3.*conj(E3)) ; i4 = R*(E4.*conj(E4)); % Quadrature;
        
        % Balanced photodetection:
        iI = i1 - i2 ; iQ  = i3 - i4;

        % Output signal:
        Out(:,3) = iI ; Out(:,4) = iQ;
    end
end