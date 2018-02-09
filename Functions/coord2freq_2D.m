function out = coord2freq_2D(T, kx, ky, w, padN)
%% coord2freq_2D.m
%
% Calculates the 2D continuous Fourier transform of T.
%
% Inputs: T      : array of 2D transfer function or signal in coordinate space 
%                  size(T) = [length(ky), length(kx), length(w)]
%         kx     : vector of spatial x frequencies
%         ky     : vector of spatial y frequencies
%         w      : vector of temporal frequencies
%         padN   : number extending the length of T to N by padding zeros
%                  (not a required input)
%
% Output: out    : array of 2D transfer function or signal in frequency space 
% 
% Example:
% >> params = loadParameters;
% >> kx = linspace(-500,500,100); ky = linspace(-500,500,100); 
% >> w = linspace(-1,1,100); 
% >> T = calcTransFuncs_fromPhi_2D(kx, ky, w, params);
% >> Tcoord = freq2coord_2D(T.T_Yphi, kx, ky, w) % gives out HRF of the 2D
%                                                 T_Yphi transfer function
% >> out = coord2freq_2D(Tcoord, kx, ky, w) % gives out 2D Fourier transform 
%                                             of the T_Yphi HRF
%
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% calculating the -1 vectors needed for the Fourier transform
kxvals = (-1).^(1:length(kx));
kyvals = (-1).^(1:length(ky));
wvals = (-1).^(1:length(w));

% converting kxvals, kyvals, and wvals into arrays
[kxM, kyM, wM] = meshgrid(kxvals, kyvals, wvals);

% performing the 2D Fourier transform 
if (nargin > 5)
    out = wM.*ifft(wM.*kxM.*fft(kxM.*kyM.*fft(kyM.*T, padN, 1), padN, 2), padN, 3);
else
    out = wM.*ifft(wM.*kxM.*fft(kxM.*kyM.*fft(kyM.*T, [], 1), [], 2), [], 3);
end;
