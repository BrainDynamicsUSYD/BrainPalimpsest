function cmap = colormap_bluetored
%% colormap_bluetored.m
%
% Makes a customized colormap from blue (most negative) to white (zero) to 
% red (most positive). 
%
% Inputs: none
%
% Output: cmap   : 80 x 3 array containing the rgb values
% 
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% total size along the first dimension of the colormap 
m = 80;

m1 = m/2;
m2 = m1 + 2;

% making the first half of the colormap (blue to white)
cmap1_r = (0:m2-1)'/(m2);
cmap1_g = cmap1_r;
cmap1_b = ones(m2,1);
cmap1 = [cmap1_r, cmap1_g, cmap1_b];

% making the second hald of the colormap (white to yellow to almost red)
cmap2_base = hot;
cmap2 = cmap2_base(end:-1:end-m1-1, :); 

% final colormap
cmap = [cmap1; cmap2];