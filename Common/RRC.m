function g = RRC(Span,SpS,Rolloff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRRC g = RRC(Span,SpS,Rolloff)                                          %
%                                                                         %
%  This function generates the impulse response of a root-raised cosine   %
% (RRC) filter with roll-off 'Rolloff', span of 'Span' symbols, and with  %
% 'SpS' samples per symbol.                                               %
%                                                                         %
% Input:                                                                  %
%   Span    = Span (in symbols) of the filter;                            %
%   SpS     = Number of samples per symbol to be considered;              %
%   Rolloff = Roll-off of the RRC filter;                                 %
%                                                                         %
% Output:                                                                 %
%   g = Impulse response of the RRC filter (column vector);               %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Configurations:
    g = zeros(Span*SpS+1,1); k = (-Span*SpS/2:Span*SpS/2)/SpS;
    i1 = find(k==0) ; i2 = find(abs(4*Rolloff*k)-1==0) ; i3 = 1:length(k);
    i3([i1 i2]) = []       ; k = k(i3);

    % Singularity in k = 0:
    if ~isempty(i1)
        g(i1) = 1-Rolloff + 4*Rolloff/pi;
    end

    % Singularity in k = 1/(4*Rolloff):
    if ~isempty(i2)
        g(i2) = Rolloff/sqrt(2)*((1+2/pi)*sin(pi/(4*Rolloff))+...
             (1-2/pi)*cos(pi/(4*Rolloff)));
    end

    % Calculating the coefficients for k ~= 0 and k ~= 1/(4*Rolloff):
    g(i3) = (sin(pi*k*(1-Rolloff)) + 4*Rolloff*k.*cos(pi*k*(1+Rolloff)))...
        ./(pi*k.*(1-(4*Rolloff*k).^2));
    
    % Normalizing the amplitude of the filter:
    g = g/max(g);
end