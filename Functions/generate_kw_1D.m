function [kx, w] = generate_kw_1D(kxsamp, wsamp, Nkx, Nw)
%% generate_kw_1D.m
%
% Generates the spatial and temporal frequency vectors.
%
% Inputs: kxsamp    : spatial sampling frequency
%         wsamp     : temporal sampling frequency
%         
%         Nkx       : number of spatial frequency points
%         Nw        : number of temporal frequency points        
%
% Output: kx        : vector of spatial frequencies
%         w         : vector of temporal frequencies        
% 
% Example:
% >> kxsamp = 0.1; wsamp = 0.1; 
% >> Nkx = 2^10; Nw = 2^10; 
% >> [kx, w] = generate_kw_1D(kxsamp, wsamp, Nkx, Nw); % gives out the frequency
%                                                  vectors
% 
% Original: Kevin Aquino, University of Sydney, 2014
%           James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% creating indices 
jveck = 0:Nkx-1;
jvec = 0:Nw-1;

% calculating the 1D spatial and temporal frequency vectors
kx = (kxsamp)*1/Nkx*(jveck - Nkx/2);
w = (wsamp)*1/Nw*(jvec - Nw/2);
