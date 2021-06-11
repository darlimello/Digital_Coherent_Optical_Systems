function [EOutput] = IQModulator(xb,EInput,ParamMZM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IQMODULATOR [EOutput] = IQModulator(xb,EInput,ParamMZM)                 %
%                                                                         %
%  This function simulates an in-phase and quadrature modulator (IQM). The%
% electrical field 'EInput' of the optical signal is modulated according  %
% to the signal 'xb', generating 'EOutput'. The Mach-Zehnder modulators   %
% (MZMs) that compose the IQM are considered to be identical and have the %
% parameters specified in 'ParamMZM'.                                     %
%                                                                         %
% Input:                                                                  %
%  xb       = (Electric) Signal to be transmitted (in one pol. orientation%
%             - column vector);                                           %
%  EInput   = Optical carrier (single pol. orientation - column vector);  %
%  ParamMZM = Struct that specifies parameters of the MZMs that compose   %
%             the IQM:                                                    %
%      -ParamMZM.Vpi    = MZM Vpi;                                        %
%      -ParamMZM.Bias   = Bias voltage;                                   %
%      -ParamMZM.MaxExc = Upper limit for the excursion of the modulation %
%                         signal*;                                        %
%      -ParamMZM.MinExc = Lower limit for the excursion of the modulation %
%                         signal*;                                        %
%       *Note: In this function, the modulation signal is scaled so that  %
%              it fits the excursion defined in 'ParamMZM.MaxExc/MinExc'; %
%                                                                         %
% Output:                                                                 %
%   EOuput = IQM output signal (column vector);                           %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Obtaining the in-phase and quadrature componentes of the electrical:
    mI = real(xb); mI = mI/max(abs(mI)); % In-phase;   
    mQ = imag(xb); mQ = mQ/max(abs(mQ)); % Quadrature;
       
    % Setting the signal excursion:
    mI = mI*(ParamMZM.MaxExc-ParamMZM.MinExc)/2; 
    mQ = mQ*(ParamMZM.MaxExc-ParamMZM.MinExc)/2;
   
    % Obtaining the signals after considering the bias:
    vI = mI+ParamMZM.Bias; vQ = mQ+ParamMZM.Bias;
        
    % Phase modulation in the in-phase and quadrature branches;
    PhiI = pi*(vI)/ParamMZM.Vpi; PhiQ = pi*(vQ)/ParamMZM.Vpi;

    % IQM output signal:
    EOutput = (0.5*cos(0.5*PhiI) + 0.5i*cos(0.5*PhiQ)).*EInput;
end