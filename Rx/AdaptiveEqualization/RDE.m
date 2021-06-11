function [w1V,w1H,w2V,w2H] = RDE(xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RDE [w1V,w1H,w2V,w2H] = RDE(xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)           %
%                                                                         %
%  This function performs the update of the filters of the MIMO butterfly %
% equalizer using the RDE algorithm.                                      %
%                                                                         %
% Input:                                                                  %
%   xV = Column vector that represents the complex samples at the MIMO    %
%        butterfly equalizer input for the vertical pol. orientation;     %
%   xH = Column vector that represents the complex samples at the MIMO    %
%        butterfly equalizer input for the horizontal pol. orientation;   %
%   y1 = Sample at the output 1 of the MIMO butterfly equalizer;          %
%   y2 = Sample at the output 2 of the MIMO butterfly equalizer;          %
%   w1V, w1H, w2V, w2H = N-coefficient FIR filters that compose the MIMO  %
%                        butterfly equalizer;                             %
%   R  = radius used as reference for coefficients calculation;           %
%   Mu = Step-size for coefficients calculation;                          %
%   Note: xV and xH must have the same length as the FIR filters w1V, w1H,%
%         w2V, w2H;                                                       %
%                                                                         %
% Output:                                                                 %
%   w1V, w1H, w2V, w2H = Updated N-coefficient FIR filters that compose   %
%                        the MIMO butterfly equalizer;                    %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Radius for output y1 and output y2:
    [~,r1] = min(abs(R-abs(y1))) ; [~,r2] = min(abs(R-abs(y2)));

    %Updating the filters:
    w1V = w1V + Mu*xV*(R(r1)^2-abs(y1).^2)*conj(y1);        
    w1H = w1H + Mu*xH*(R(r1)^2-abs(y1).^2)*conj(y1);   
    w2V = w2V + Mu*xV*(R(r2)^2-abs(y2).^2)*conj(y2);  
    w2H = w2H + Mu*xH*(R(r2)^2-abs(y2).^2)*conj(y2);    
end