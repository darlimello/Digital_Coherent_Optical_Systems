function [xb] = PulseShaping(x,SpS,ParamFilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PULSESHAPING xb = PulseShaping(x,SpS,ParamFilter)                       %
%                                                                         %
%  This function performs pulse shaping of a sequence of symbols 'x' using%
% a root-raised cosine (RRC) filter. The filter parameters are defined in %
% 'ParamFilter'. The sequence of symbols is first upsampled to 'SpS'      %
% samples per symbol and, then, applied to the RRC filter.                %
%                                                                         %
% Input:                                                                  %
%   x           = Sequence of symbols to be transmitted (one pol.         %
%                 orientation - column vector);                           %
%   SpS         = Number of samples per symbol to be considered during    %
%                 pulse shaping;                                          %
%   ParamFilter = Struct that specifies parameters of the RRC filter:     %
%          - ParamFilter.Rolloff: Roll-off factor (between 0 and 1);      %
%          - ParamFilter.Span: Filter span (in symbols);                  %
%                                                                         %
% Output:                                                                 %
%   xb = Signal after pulse shaping, normalized to unitary power (column  %
%        vector);                                                         %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Obtaining the filter transfer function:   
    g = RRC(ParamFilter.Span,SpS,ParamFilter.Rolloff);

    % Upsampling the symbols for pulse shaping:
    xUpsamp = upsample(x(:,1),SpS);
    
    % Filtering the upsampled symbols:
    xb(:,1) = conv(xUpsamp,g,'same');
    
    % Normalizing the signal to unitary power:
    xb(:,1) = xb(:,1)/sqrt(mean(abs(xb(:,1)).^2));    
end