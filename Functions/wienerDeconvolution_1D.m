function deconvolved = wienerDeconvolution_1D(signal, trans_func, ...
                                           kx, w, noise)
%% wienerDeconvolution_1D.m
%
% Performs the 1D Wiener Deconvolution of the response, corresponding to 
% trans_func, from signal.
%
% Inputs: signal        : array of 1D signal (x,t)
%                         size(signal) = [length(kx), length(w)]                          
%         trans_func    : array of 1D transfer function corresponding to the 
%                         response to be deconvolved
%                         size(trans_func) = [length(kx), length(w)]   
%         kx            : vector of spatial frequencies
%         w             : vector of temporal frequencies         
%         noise         : if length(noise) > 1
%                            noise(1) is the spatial cut-off frequency
%                            noise(2) is the temporal cut-off frequency
%                         if length(noise) == 1
%                            noise is the constant sigma 
%
% Output: deconvolved   : structure containing the deconvolution result.
%                         Possible fields are coord (for coordinate space)
%                         and freq (for frequency space).
% 
% Example:
% >> params = loadParameters;
% >> kx = linspace(-500,500,100); w = linspace(-1,1,100); 
% >> T = calcTransFuncs_fromPhi_1D(kx, w, params);
% >> load BOLD_signal.mat   % assuming the data is stored in this mat file
% >> noise = 0.1;           % assuming a constant NSR
% >> deconvolved = wienerDeconvolution_1D(BOLD_signal, T.T_Yphi, kx, w, noise) 
% >> deconvolved.coord      % gives out the deconvolved phi in coordinate space
% 
% Kevin Aquino, University of Sydney, 2014
% James Pang, University of Sydney, 2016

%%

% converting the signal to frequency space
signal_freq =  coord2freq_1D(signal, kx, w);

% calculating the noise-to-signal ratio (NSR)
if (length(noise)>1)
    [~, kxCloseInd] = min(abs(kx - noise(1)*2*pi));
    [~, wCloseInd] = min(abs(w - noise(2)*2*pi));
    
    NSR = trans_func(kxCloseInd, wCloseInd);
else
    NSR = noise;
end;

% calculating the Wiener filter
filter = (1./trans_func).*(abs(trans_func).^2./(abs(trans_func).^2 + abs(NSR)^2));

% calculating the deconvolved response in frequency space
deconvolved.freq = filter.*signal_freq;

% calculating the deconvolved response in coordinate space
deconvolved.coord = freq2coord_1D(deconvolved.freq, kx, w);
