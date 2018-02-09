function [overlay_quantity] = goingBackToTriangulation(hemisphere, scanNo, ...
                            grid_quantity, interpolant, xcoords, ycoords)
%% goingBackToTriangulation.m
%
% Converts the flattened gridded BOLD quantity to flattened triangulated
% form compatible as an overlay for freesurfer.
%
% Inputs: hemisphere        : string of hemisphere
%                             Possible fields are lh for left hemisphere 
%                             and rh for right hemisphere.
%         scanNo            : scan number to be processed
%                             Possible fields are 
%                             1 for the expanding ring data
%                             2, 3, ..., 11 for the expanding arc data.
%         grid_quantiy      : array of gridded quantity 
%         F                 : interpolant object
%         xcoords           : vector of x coordinates
%         ycoords           : vector of y coordinates
%
% Output: overlay_quantity  : converted triangulated quantity 
% 
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%% Preparing the vertices, faces, and overlay

% reading vertices of flat surface
flat = read_patch(['Data/ExpandingRingAndExpandingArc/FreesurferFiles/',hemisphere,...
                       '.occip.flat.patch.3d']);

% reading the vertex coordinates and face lists of white matter surface
[~, fac] = read_surf(['Data/ExpandingRingAndExpandingArc/FreesurferFiles/',hemisphere,...
                       '.white']);

% finding the face index of flat patch corresponding to white surface
[h, h2] = ismember(fac, flat.ind);
sum3 = sum(h, 2);
fac2 = h2(sum3==3, :);

% constructing the vertex coordinates and faces of flat patch in the actual
% brain for visualization
flat_struct.Vertices = [flat.x; flat.y; flat.z].';

% loading an overlay: BOLD data
overlay_BOLD = MRIread(['Data/ExpandingRingAndExpandingArc/OriginalSurfaces/', ...
                         hemisphere,'.Scan',num2str(scanNo),'.mgz']);

%% Getting the triangulated overlays

overlay_quantity = overlay_BOLD;
overlay_quantity.vol = NaN(size(overlay_BOLD.vol));

% check if xcoords or ycoords are in m or mm
if max(xcoords) < 1 
    Kx = dsearchn((xcoords*1e3).', interpolant.Points(:,1));
else
    Kx = dsearchn(xcoords.', interpolant.Points(:,1));
end
if max(ycoords) < 1 
    Ky = dsearchn((ycoords*1e3).', interpolant.Points(:,2));
else
    Ky = dsearchn(ycoords.', interpolant.Points(:,2));
end

for i = 1:size(grid_quantity, 3)
    quantity_vol = zeros(1, size(interpolant.Values, 1)); 
    for j = 1:size(interpolant.Values, 1)
    	quantity_vol(1, j) = grid_quantity(Ky(j), Kx(j), i);
    end
    
    overlay_quantity.vol(1, flat.ind+1, 1, i) = quantity_vol;
end

overlay_quantity.vol = overlay_quantity.vol(:,:,:,1:size(grid_quantity, 3));
