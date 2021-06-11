function [Out] = InsertSkew(In,SpS,Rs,NPol,ParamSkew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSERTSKEW [Out] = InsertSkew(In,SpS,Rs,NPol,ParamSkew)                 %
%                                                                         %
%   This function inserts a temporal misalignment (skew) between in-phase %
% and quadrature components of the signal 'In'. Signal 'In' must be the   %
% signal obtained at the output of the optical front-end, right before    %
% the 'ADC'. In this function, a temporal delay is specified for each     %
% component in 'ParamSkew'. The temporal misalignment between the         %
% components is then applied assuming the minimum temporal delay as       %
% reference. E.g.: A skew of 10 ps between the in-phase and quadrature    %
% components of the signal in V pol. orientation can be specified using   %
% 'ParamSkew.TauIV = 5e-12' and 'ParamSkew.TauQV = -5e-12'. For           %
% transmission with pol. multiplexing, a unique reference is used for the %
% four components (in-phase and quadrature components of each pol.        %
% orientation).                                                           %
%                                                                         %
% Input:                                                                  %
%   In        = Input signals. For transmission in single pol. orientation%
%               'In' must be a matrix with two column vectors. The 1st and%
%               2nd columns must have the in-phase and quadrature         %
%               components of the signal (e.g., at the output of the ADC).%
%               For transmission with pol. multiplexing, 'In' must be a   %
%               matrix with four column vectors. The 1st and the 2nd      %
%               columns must have the in-phase and quadrature components  %
%               of the signal received in the V pol. orientation,         %
%               respectively. The 3rd and 4th columns have the in-phase   %
%               and quadrature components of the signal received in the H %
%               pol. orientation, respectively;                           %
%   SpS       = Number of samples per symbol in the input signal 'In';    %
%   Rs        = Symbol rate in Symbols/second;                            %
%   NPol      = Number of polarization orientations;                      %
%   ParamSkew = Struct that specifies the temporal delay (in seconds) for %
%               each component of the input signal, i.e,                  %
%                   ParamSkew.TauIV and ParamSkew.TauQV;                  %
%                   ParamSkew.TauIH and ParamSkew.TauQH (if NPol == 2);   %
%   *Note: 'In' is the signal obtained at the output of the optical front-%
%           end, right before the 'ADC'.                                  %
%                                                                         %
% Output:                                                                 %
%   Out = Signal produced after skew insertion. 'Out' is organized in     %
%         columns in the same way as 'In'. E.g., 'Out' is the signal at   %
%         the input of the ADC;                                           %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input signals:
    iIV = In(:,1) ; iQV = In(:,2);
    if NPol == 2
        iIH = In(:,3) ; iQH = In(:,4);
    end
    
    % Calculating the interpolation factor given the timing skew
    % (TauIV/H and TauQV/H, for in-phase and quadrature components of 
    % V pol. and H pol. orientations, respectively):
    Ts = 1/(SpS*Rs) ; Skew = [ParamSkew.TauIV/Ts ParamSkew.TauQV/Ts]; 
    if NPol == 2
        Skew = [Skew ParamSkew.TauIH/Ts ParamSkew.TauQH/Ts];       
    end

    % Using the min skew as reference:
    Skew = Skew-min(Skew);

    % Inserting skew in the samples:
    Len = length(iIV);
    iIV = interp1(0:Len-1,iIV,Skew(1):Len-1,'spline','extrap').';
    iQV = interp1(0:Len-1,iQV,Skew(2):Len-1,'spline','extrap').';
    if NPol == 2
        iIH = interp1(0:Len-1,iIH,Skew(3):Len-1,'spline','extrap').';
        iQH = interp1(0:Len-1,iQH,Skew(4):Len-1,'spline','extrap').';
        
        % Output signals with the same length:
        MinLength = min([length(iIV) length(iQV) length(iIH) length(iQH)]);
        Out(:,1) = iIV(1:MinLength) ; Out(:,2) = iQV(1:MinLength);
        Out(:,3) = iIH(1:MinLength) ; Out(:,4) = iQH(1:MinLength);
    else
        % Output signals with the same length:
        MinLength = min([length(iIV) length(iQV)]);
        Out(:,1)  = iIV(1:MinLength) ; Out(:,2) = iQV(1:MinLength);
    end
end