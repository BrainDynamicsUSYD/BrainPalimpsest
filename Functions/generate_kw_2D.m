function [kx, ky, w] = generate_kw_2D(kxsamp, kysamp, wsamp, Nkx, Nky, Nw)
%% generate_kw_2D.m
%
% Generates the spatial and temporal frequency vectors.
%
% Inputs: kxsamp    : spatial x sampling frequency
%         kysamp    : spatial y sampling frequency
%         wsamp     : temporal sampling frequency
%         
%         Nkx       : number of spatial x frequency points
%         Nky       : number of spatial y frequency points
%         Nw        : number of temporal frequency points        
%
% Output: kx        : vector of spatial x frequencies
%         ky        : vector of spatial y frequencies
%         w         : vector of temporal frequencies        
% 
% Example:
% >> kxsamp = 0.1; kysamp = 0.1; wsamp = 0.1; 
% >> Nkx = 2^10; Nky = 2^10; Nw = 2^10; 
% >> [kx, ky, w] = generate_kw_2D(kxsamp, kysamp, wsamp, Nkx, Nky, Nw); 
                                          % gives out the frequency vectors
% 
% James Pang, University of Sydney, 2016

%%

% creating indices 
jveckx = 0:Nkx-1;
jvecky = 0:Nky-1;
jvec = 0:Nw-1;

% calculating the spatial and temporal frequency vectors
kx = (kxsamp)*1/Nkx*(jveckx - Nkx/2);
ky = (kysamp)*1/Nky*(jvecky - Nky/2);
w = (wsamp)*1/Nw*(jvec - Nw/2);
