function out = freq2coord_1D(T, kx, w, padN)
%% freq2coord_1D.m
%
% Calculates the continuous 1D inverse Fourier transform of T.
%
% Inputs: T      : array of 1D transfer function or signal in frequency space 
%                  size(T) = [length(kx), length(w)]
%         kx     : vector of spatial frequencies
%         w      : vector of temporal frequencies      
%         padN   : number extending the length of T to N by padding zeros
%                  (not a required input)
%
% Output: out    : array of transfer function or signal in coordinate space 
% 
% Example:
% >> params = loadParameters;
% >> kx = linspace(-500,500,100); w = linspace(-1,1,100); 
% >> T = calcTransFuncs_fromPhi_1D(kx, w, params);
% >> out = freq2coord_1D(T.T_Yphi, kx, w) % gives out inverse Fourier transform 
%                                           of the T_Yphi transfer function
%
% Kevin Aquino, University of Sydney, 2014
% Original: Kevin Aquino, University of Sydney, 2014
%           James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% calculating the -1 vectors needed for the Fourier transform
kxvals = (-1).^(1:length(kx));
wvals = (-1).^(1:length(w));

% converting kxvals and wvals into arrays
kxM = repmat(kxvals.', 1, length(w));
wM = repmat(wvals, length(kx), 1);

% performing the Fourier transform 
if (nargin > 5)
    out = (kxM.*ifft(kxM.*wM.*fft(wM.*T, padN, 2), padN, 1));
else
    out = (kxM.*ifft(kxM.*wM.*fft(wM.*T, [], 2), [], 1));
end;
