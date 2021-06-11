function [y] = AdaptiveEqualizer(x,SpS,ParamDE)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAPTIVEEQUALIZER [y] = AdaptiveEqualizer(x,SpS,ParamDE)                %
%                                                                         %
%   This function performs equalization in a dual polarization signal 'y' %
% using the CMA or the RDE algorithms. The supported modulation formats   %
% are QPSK and 16-QAM. The CMA can be used for pre-convergence of the RDE %
% algorithm.                                                              %
%                                                                         %
% Input:                                                                  %
%  x       = Input signal for two pol. orientation (matrix with two       %
%            column vectors, where each column corresponds to the signal  %
%            of one pol. orientation). 'x' must be normalized to unitary  %
%            power and obtained at 2 Sa/s;                                %
%  SpS     = Number of samples per symbol at the input signal;            %
%  ParamDE =Struct that specifies parameters for the adaptive equalization%
%      - ParamDE.Eq: Defines the algorithm to be used:                    %
%                    'CMA' - CMA only;                                    %
%                    'RDE' - RDE only;                                    %
%                    'CMA+RDE' - CMA is used to initialize the RDE;       %
%      - ParamDE.NTaps: Number of taps for the filters in the butterfly   %
%                       configuration;                                    %
%      - ParamDE.Mu: Step-size for coefficients calculation;              %
%      - ParamDE.SingleSpike: 'true': Single spike initialization;        %
%                             'false': All taps are initialized with zeros%
%      - ParamDE.N1: Number of coefficient calculations to perform prior  %
%                    to proper initialization of the filters w2H and w2V; %
%            ('ParamDE.N1' is related to the single spike initialization) %
%      - ParamDE.N2: Number of coefficient calculations to perform prior  %
%                    to swicth from CMA to RDE (Note that N2 must only be %
%                    defined if CMA is used for RDE initialization);      %
%      - ParamDE.NOut: Number of samples to discard after equalization;   %
%                                                                         %
% Output:                                                                 %
%   y = Output signal (at 1 Sa/s) after adaptive equalization;            %
%                                                                         %
% This function is part of the book Digital Coherent Optical Systems;     %
% Darli A. A. Mello and Fabio A. Barbosa;                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Equalization algorithms:
    CMAFlag = false ; RDEFlag = false ; CMAtoRDE = false ; CMAInit = false;
    if strcmp(ParamDE.Eq,'CMA')
        CMAFlag  = true;
    elseif strcmp(ParamDE.Eq,'RDE')
        RDEFlag  = true;
    elseif strcmp(ParamDE.Eq,'CMA+RDE')
        CMAFlag = true ; CMAtoRDE = true ; CMAInit = true;
    end        
    
    % Radii for constellations with unitary power:
    if CMAFlag
        % Radius for the CMA (E{|x|^4}/E{|x|^2}):
        if ~CMAtoRDE
            R_CMA = 1;
        else
            R_CMA = 1.32;
        end
    end
    if CMAtoRDE || RDEFlag       
        R_RDE = [1/sqrt(5) 1 3/sqrt(5)];
    end
    
    % Important parameters:
    Mu    = ParamDE.Mu    ; SingleSpike = ParamDE.SingleSpike; 
    NTaps = ParamDE.NTaps ; N1 = ParamDE.N1 ; NOut = ParamDE.NOut;
    if isfield(ParamDE,'N2') && CMAInit
        N2 = ParamDE.N2;
    end
            
    % Input blocks:
    x  = [x(end-floor(NTaps/2)+1:end,:) ; x ; x(1:floor(NTaps/2),:)];
    xV = convmtx(x(:,1).',NTaps)     ; xH = convmtx(x(:,2).',NTaps);
    xV = xV(:,NTaps:SpS:end-NTaps+1) ; xH = xH(:,NTaps:SpS:end-NTaps+1);
    
    % Output length:
    OutLength = floor((size(x,1)-NTaps+1)/2) ; clearvars x

    % Initializing the outputs
    y1 = zeros(OutLength,1) ;  y2 = zeros(OutLength,1);
    
    % Initial filter coefficients:
    w1V = zeros(NTaps,1); w1H = zeros(NTaps,1); w2V = zeros(NTaps,1);
    w2H = zeros(NTaps,1);
    
    % If single spike initialization:
    if SingleSpike
       w1V(floor(NTaps/2)+1) = 1;
    end
    
    for i = 1:OutLength    
        % Calculating the outputs:
        y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
        y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);
        
        % Updating the filter coefficients:
        if CMAFlag            
            % Constant modulus algorithm:
            [w1V,w1H,w2V,w2H] = CMA(xV(:,i),xH(:,i),y1(i),y2(i),w1V,...
                w1H,w2V,w2H,R_CMA,Mu);     
            
            % Switching from 
            if CMAtoRDE
                if i == N2
                    CMAFlag = false ; RDEFlag = true;
                end
            end             
        elseif RDEFlag            
            % Radius-directed equalization:
            [w1V,w1H,w2V,w2H] = RDE(xV(:,i),xH(:,i),y1(i),y2(i),w1V,...
                w1H,w2V,w2H,R_RDE,Mu);            
        end
    
        % Reinitialization of the filter coefficients:
        if and(i == N1,SingleSpike)
            w2H = conj(w1V(end:-1:1,1)) ; w2V = -conj(w1H(end:-1:1,1));            
        end               
    end
    
    % Output samples:
    y = [y1 y2] ; y = y(1+NOut:end,:);
end