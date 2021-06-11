function [r] = ADC(i,SpSIn,NPol,ParamFilter,ParamADC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC [r] = ADC(i,SpSIn,NPol,ParamFilter,ParamADC)                        %
%                                                                         %
%  This function downsamples the signal 'i' so that the output signal 'r' %
% has 'ParamADC.SpS' samples per symbol. Before downsampling, signal 'i'  %
% is filtered by a filter defined in 'ParamFilter'.                       %
%                                                                         %
% Input:                                                                  %
%   i = Signal to be sampled. For transmission in single pol. orientation,%
%       'i' must be a matrix with two column vectors. The 1st column must %
%       have the in-phase component of the signal, while the 2nd column   %
%       must have the quadrature component. For transmission with pol.    %
%       multiplexing, 'i' must be a matrix with four column vectors. The  %
%       1st and 2nd columns must have the in-phase and quadrature         %
%       components of the signal in the V pol. orientation, respectively. %
%       The 3rd and 4th columns must have the in-phase and quadrature     %
%       components of the signal in the H pol. orientation, respectively; %
%   SpSIn       = Number of samples per symbol at the input signal 'i';   %
%   NPol        = Number of polarization orientations used;               %
%   ParamFilter = Struct that specifies parameters of the filter used     %
%                 before sampling 'i', e.g. 'ParamFilter.Type', which can %
%                 be 'RRC','SuperGaussian', or 'NoFiter'.                 %
%     -If ParamFilter.Type = 'RRC', additional required parameters are:   %
%         *ParamFilter.Rolloff: Roll-off of the RRC filter;               %
%         *ParamFilter.Span: Span (in symbols) of the RRC filter          %
%     -If ParamFilter.Type = 'SuperGaussian', additional required         %
%      parameters are:                                                    %
%         *ParamFilter.Order: Order of the SuperGaussian filter           %
%         *ParamFilter.Bw: (Baseband) Filter bandwidth normalized to the  %
%                          symbol rate.                                   %
%     -If ParamFilter.Type = 'NoFilter', the signals are not filtered     %
%         prior to downsampling.                                          %
%   ParamADC    = Struct that specifies parameters of the ADC:            %
%     - ParamADC.SpS = Number of samples per symbol required at the output%
%       signal 'r';                                                       %
%     - ParamADC.FreqError = Deviation in parts per million (ppm) of the  %
%       sampling rate. The default is 0 ppm (no sampling errors);         %
%     - ParamADC.PhaseError = Constant lag with the optimal sampling      %
%       instant normalized by the symbol period Ts. Its interval is       %
%       limited to [-Ts/2,Ts/2]. The default value is 0;                  %
%                                                                         %
% Output:                                                                 %
%   r = Signal after filtering and downsampling to 'ParamADC.SpS' samples %
%       per symbol. 'r' is organized in columns in the same way as 'i';   %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~strcmp(ParamFilter.Type,'NoFilter')
        switch ParamFilter.Type
            case 'RRC'
                % Obtaining the filter transfer function:
                g = RRC(ParamFilter.Span,SpSIn,ParamFilter.Rolloff);

                % Filtering the signals:
                i(:,1)=conv(i(:,1),g,'same'); i(:,2)=conv(i(:,2),g,'same');
                if NPol == 2
                i(:,3)=conv(i(:,3),g,'same'); i(:,4)=conv(i(:,4),g,'same');
                end 
            case 'SuperGaussian'
                % Frequency vector
                f = fftshift(-0.5:1/length(i(:,1)):0.5-1/length(i(:,1)))...
                    .'*SpSIn;

                % Super-Gaussian filter in the frequency domain:
                G = exp(-log(sqrt(2))*(f/ParamFilter.Bw)...
                    .^(2*ParamFilter.Order));

                % Filtering the signals:
                i(:,1)=ifft(fft(i(:,1)).*G) ; i(:,2)=ifft(fft(i(:,2)).*G);
                if NPol == 2
                  % Filtering the signals:
                  i(:,3)=ifft(fft(i(:,3)).*G); i(:,4)=ifft(fft(i(:,4)).*G);
                end
            otherwise
                error('Filter type not supported');            
        end
    end
    
    % Normalizing the signals to unitary energy:
    i(:,1) = i(:,1)/sqrt(mean(abs(i(:,1)).^2));
    i(:,2) = i(:,2)/sqrt(mean(abs(i(:,2)).^2));    
    if NPol == 2
        i(:,3) = i(:,3)/sqrt(mean(abs(i(:,3)).^2));
        i(:,4) = i(:,4)/sqrt(mean(abs(i(:,4)).^2));
    end
   
    % Downsampling the signals to 2 samples per symbol:
    PhaseError = 0; % Default value;
    if isfield(ParamADC,'PhaseError')
        PhaseError = ParamADC.PhaseError*SpSIn;
    end    
    Len = length(i(:,1)) ; Pos = PhaseError+(1:(SpSIn/ParamADC.SpS):Len)';
    r(:,1) = interp1(1:Len,i(:,1),Pos,'spline','extrap');
    r(:,2) = interp1(1:Len,i(:,2),Pos,'spline','extrap');
    if NPol == 2
        % Downsampling the signals to 2 samples per symbol:
        r(:,3) = interp1(1:Len,i(:,3),Pos,'spline','extrap');
        r(:,4) = interp1(1:Len,i(:,4),Pos,'spline','extrap');
    end    
    
   if isfield(ParamADC,'FreqError')
       Len = length(r(:,1)) ; Pos = (1:1-ParamADC.FreqError:Len)';
       rAux(:,1) = interp1(1:Len,r(:,1),Pos,'spline','extrap');
       rAux(:,2) = interp1(1:Len,r(:,2),Pos,'spline','extrap');
       if NPol == 2
           rAux(:,3) = interp1(1:Len,r(:,3),Pos,'spline','extrap');
           rAux(:,4) = interp1(1:Len,r(:,4),Pos,'spline','extrap'); 
       end
       r = rAux;
   end   
end