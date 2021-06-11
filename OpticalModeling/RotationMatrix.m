function [Out] = RotationMatrix(In,TDegree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATIONMATRIX [Out] = RotationMatrix(In,TDegree)                       %
%                                                                         %
%   This function mixes the signals of a polarization-multiplexed         %
% transmission applying a rotation matrix                                 %
%                 R = [cos(Theta) -sin(Theta)                             %
%                      sin(Theta)  cos(Theta)],                           % 
% where Theta is the rotation angle in radians;                           %
%                                                                         %
% Input:                                                                  %
%   In      = Input signal. 'In' must be a matrix with two columns, where %
%             each column has the signal of a polarization orientation (V %
%             and H pol. orientations);                                   %
%   TDegree = Rotation angle in degrees to be considered in the rotation  %
%             matrix;                                                     %
%                                                                         %
% Output:                                                                 %
%   Out = Output signals after mixing using the rotation matrix 'R'. 'Out'%
%         is organized in columns in the same way as 'In';                %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Rotation angle in radians:
    Theta = TDegree*pi/180;

    % Rotation matrix:
    R = [cos(Theta) -sin(Theta); sin(Theta)  cos(Theta)];

    % Applying the rotation to the signals:
    Out(:,1) = R(1,1)*In(:,1) + R(1,2)*In(:,2); 
    Out(:,2) = R(2,1)*In(:,1) + R(2,2)*In(:,2); 
end