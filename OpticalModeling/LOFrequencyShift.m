function [EOut] = LOFrequencyShift(EIn,Delta_f,SpS,Rs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOInsertion [EOut] = LOFrequencyShift(EIn,Delta_f,SpS,Rs)               %
%                                                                         %
%  This function shifts the carrier frequency of the local oscillator (LO)%
% laser by '-Delta_f'. To do so, a constant phase offset between          %
% neighboring samples calculated as                                       %
%                   DeltaTheta_f = -2*pi*Delta_f*T                        %
% where 'T' is the sampling period (in seconds) and 'Delta_f' is the      %
% carrier frequency shift (in Hz), is applied to the signal generated by  %
% the LO laser, 'EIn'. This function can be used to insert a carrier      %
% frequency offset in the coherent detected signals. Assuming that the    %
% carrier frequency of the transmitter and LO lasers are 'X Hz', shifting %
% the carrier frequency of the LO laser by '-Delta_f' produces a carrier  %
% frequency offset of 'Delta_f' after coherent detection.                 % 
%                                                                         %
% Input:                                                                  %
%   EIn     = Signal produced by the LO laser (e.g. produced by  the      %
%             function 'Laser'). For transmission in single pol.          %
%             orientation, 'EIn' must be a column vector. For transmission%
%             with pol. multiplexing, 'EIn' must be a matrix with two     %
%             column-oriented vectors, where each column vector represents%
%             the signal of a pol. orientation;                           %
%   Delta_f = Carrier frequency shift (in Hz) to be applied in the signal %
%             'EIn';                                                      % 
%   SpS     = Number of samples per symbol at the oversampled signal      %
%             (transmitted signal / received signal before sampling);     %
%   Rs      = Symbol rate in symbols/second;                              %
%                                                                         %
% Output:                                                                 % 
%   EOut = Signal generated after applying the frequency shift in 'Ein';  %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Period between samples at the (oversampled) transmitted signal:
    T = 1/(SpS*Rs);

    % Phase offset between consecutive samples due to 'Delta_f':
    DeltaTheta_f = -2*pi*Delta_f*T; 

    % Inserting the effects of frequency offset into signal 'EIn':
    k    = repmat((0:size(EIn,1)-1).',1,size(EIn,2));
    EOut = EIn.*exp(1i*k*DeltaTheta_f);
end