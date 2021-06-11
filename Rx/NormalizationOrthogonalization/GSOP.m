function [rOut] = GSOP(rIn,NPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GSOP [rOut] = GSOP(rIn,NPol)                                            %
%                                                                         %
%   This function performs Gram-Schmidt orthogonalization in signal 'rIn'.%
%                                                                         %
% Input:                                                                  %
%   rIn = Input signal in which orthogonalization will be performed.      %
%         For transmission in single pol. orientation, 'In' must be a     %
%         matrix with two columns. The 1st and 2nd columns must have the  %
%         in-phase and quadrature components of the signal, respectively. %
%         For transmission with pol. multiplexing, 'In' must be a matrix  %
%         with four columns. The 1st and the 2nd columns must have the    %
%         in-phase and quadrature components of the signal in V pol.      %
%         orientation, while the 3rd and 4th columns, the in-phase and    %
%         quadrature components of the signal in H pol. orientation,      %
%         respectively. 'rIn' is, e.g., the signal at the output of the   %
%         'Deskew' stage;                                                 %
%   NPol = Number of polarization orientations used;                      % 
%                                                                         %
% Output:                                                                 %
%   rOut = Signal after Gram-Schmidt orthogonalization. 'rOut' is arranged%
%          in columns in the same way as 'rIn'. Each signal is normalized %
%          to unitary power;                                              %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Taking the in-phase component as reference:
    rIOrt = rIn(:,1)/sqrt(mean(rIn(:,1).^2));

    % Orthogonalization:
    rQInt = rIn(:,2)-mean(rIn(:,1).*rIn(:,2))*rIn(:,1)/mean(rIn(:,1).^2);
    rQOrt = rQInt/sqrt(mean(rQInt.^2));

    % Complex output signal:
    rOut = [rIOrt rQOrt];  
    
    if NPol == 2
      % Taking the in-phase component as reference:
      rIOrt = rIn(:,3)/sqrt(mean(rIn(:,3).^2));

      % Orthogonalization:
      rQInt = rIn(:,4)-mean(rIn(:,3).*rIn(:,4))*rIn(:,3)/mean(rIn(:,3).^2);
      rQOrt = rQInt/sqrt(mean(rQInt.^2));

      % Complex output signal:
      rOut = [rOut rIOrt rQOrt];
    end  
end