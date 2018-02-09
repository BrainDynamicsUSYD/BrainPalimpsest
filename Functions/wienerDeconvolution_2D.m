function deconvolved = wienerDeconvolution_2D(signal, trans_func, ...
                                           kx, ky, w, noise)
%% wienerDeconvolution_2D.m
%
% Performs the 2D Wiener Deconvolution of the response, corresponding to 
% trans_func, from signal.
%
% Inputs: signal        : array of 2D signal (x,y,t)
%                         size(signal) = [length(ky), length(kx), length(w)]                          
%         trans_func    : array of 2D transfer function corresponding to the 
%                         response to be deconvolved
%                         size(trans_func) = [length(ky), length(kx), length(w)]   
%         kx            : vector of spatial x frequencies
%         ky            : vector of spatial y frequencies
%         w             : vector of temporal frequencies         
%         noise         : if length(noise) > 1
%                            noise(1) is the spatial x cut-off frequency
%                            noise(2) is the spatial y cut-off frequency
%                            noise(3) is the temporal cut-off frequency
%                         if length(noise) == 1
%                            noise is the constant sigma 
%
% Output: deconvolved   : structure containing the deconvolution result.
%                         Possible fields are coord (for coordinate space)
%                         and freq (for frequency space).
% 
% Example:
% >> params = loadParameters;
% >> kx = linspace(-500,500,100); ky = linspace(-500,500,100); w = linspace(-1,1,100); 
% >> T = calcTransFuncs_fromPhi_2D(kx, ky, w, params);
% >> load BOLD_signal.mat   % assuming the data is stored in this mat file
% >> noise = 0.1;           % assuming a constant NSR
% >> deconvolved = wienerDeconvolution_2D(BOLD_signal, T.T_Yphi, kx, ky, w, noise) 
% >> deconvolved.coord      % gives out the deconvolved phi in coordinate space
% 
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% converting the signal to frequency space
signal_freq =  coord2freq_2D(signal, kx, ky, w);

% calculating the noise-to-signal ratio (NSR)
if (length(noise)>1)
    [~, kxCloseInd] = min(abs(kx - noise(1)*2*pi));
    [~, kyCloseInd] = min(abs(ky - noise(2)*2*pi));
    [~, wCloseInd] = min(abs(w - noise(3)*2*pi));
    
    NSR = trans_func(kyCloseInd, kxCloseInd, wCloseInd);
else
    NSR = noise;
end;

% calculating the Wiener filter
filter = (1./trans_func).*(abs(trans_func).^2./(abs(trans_func).^2 + abs(NSR)^2));

% calculating the deconvolved response in frequency space
deconvolved.freq = filter.*signal_freq;

% calculating the deconvolved response in coordinate space
deconvolved.coord = freq2coord_2D(deconvolved.freq, kx, ky, w);
