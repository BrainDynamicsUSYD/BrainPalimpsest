function [grid_Benson, F_Benson] = makingGriddedFlat_Benson(hemisphere, ...
                                                   resolution, saveMATfile)
%% makingGriddedFlat_Benson.m
%
% Converts the flattened triangulated Benson template to flattened gridded 
% (matrix) form.
%
% Inputs: hemisphere    : string of hemisphere
%                         Possible fields are lh for left hemisphere and 
%                         rh for right hemisphere.
%         resolution    : x and y spatial resolution in mm
%         saveMatfile   : 1 or 0
%                         Choose 1 if you want to save the mat file of the 
%                         converted gridded Benson map, 0 otherwise.
%
% Output: grid_Benson   : matrix of converted gridded Benson map
%         F_Benson      : interpolant object
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
flat_struct.Faces = fac2;

% loading an overlay: Benson map
overlay_Benson = MRIread(['Data/ExpandingRingAndExpandingArc/FreesurferFiles/',hemisphere,...
                       '.bensonmap.mgz']);

%% Interpolating Benson map

xcoords = flat_struct.Vertices(:, 1);
ycoords = flat_struct.Vertices(:, 2);

x_lim = floor(max(xcoords));
y_lim = floor(max(ycoords));
new_xcoords = -x_lim:resolution:x_lim;
new_ycoords = -y_lim:resolution:y_lim;
[mesh_xcoords, mesh_ycoords] = meshgrid(new_xcoords, new_ycoords);

grid_Benson = zeros(length(new_ycoords), length(new_xcoords), size(overlay_Benson.vol, 4));

for i = 1:size(overlay_Benson.vol, 4)
    Benson = squeeze(overlay_Benson.vol(1, flat.ind+1, 1, i));
    if i==1
        F_Benson = scatteredInterpolant(xcoords, ycoords, Benson');
        F_Benson.Method = 'linear';
        F_Benson.ExtrapolationMethod = 'none';
    else
        F_Benson.Values = Benson';
    end
            
    grid_Benson(:,:,i) = F_Benson(mesh_xcoords, mesh_ycoords);
end

%% Saving grid_Benson and F_Benson in a mat file

if saveMATfile
    filename = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/',hemisphere,...
                '.Benson_resolution=',num2str(resolution),'.mat'];
    save(filename, 'grid_Benson', 'F_Benson')
end
 