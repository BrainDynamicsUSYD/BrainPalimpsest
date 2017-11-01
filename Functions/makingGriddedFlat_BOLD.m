function [grid_BOLD, F] = makingGriddedFlat_BOLD(hemisphere, ...
                                            scanNo, resolution, saveMATfile)
%% makingGriddedFlat_BOLD.m
%
% Converts the flattened triangulated BOLD data to flattened gridded (matrix) 
% form.
%
% Inputs: hemisphere    : string of hemisphere
%                         Possible fields are lh for left hemisphere and 
%                         rh for right hemisphere.
%         scanNo        : scan number to be processed
%                         Possible fields are 
%                         1 for the expanding ring data
%                         2, 3, ..., 11 for the expanding arc data.
%         resolution    : x and y spatial resolution in mm
%         saveMatfile   : 1 or 0
%                         Choose 1 if you want to save the mat file of the 
%                         converted gridded BOLD data, 0 otherwise.
%
% Output: grid_BOLD     : matrix of converted gridded BOLD data
%         F             : interpolant object
% 
% James Pang, University of Sydney, 2017

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

% loading an overlay: BOLD data
overlay_BOLD = MRIread(['Data/ExpandingRingAndExpandingArc/OriginalSurfaces/', ...
                         hemisphere,'.Scan',num2str(scanNo),'.mgz']);

%% Interpolating BOLD data

xcoords = flat_struct.Vertices(:, 1);
ycoords = flat_struct.Vertices(:, 2);

x_lim = floor(max(xcoords));
y_lim = floor(max(ycoords));
new_xcoords = -x_lim:resolution:x_lim;
new_ycoords = -y_lim:resolution:y_lim;
[mesh_xcoords, mesh_ycoords] = meshgrid(new_xcoords, new_ycoords);

grid_BOLD = zeros(length(new_ycoords), length(new_xcoords), size(overlay_BOLD.vol, 4));

for i = 1:size(overlay_BOLD.vol, 4)
    BOLD = squeeze(overlay_BOLD.vol(1, flat.ind+1, 1, i));
    if i==1
        F = scatteredInterpolant(xcoords, ycoords, BOLD');
        F.Method = 'linear';
        F.ExtrapolationMethod = 'none';
    else
        F.Values = BOLD';
    end
            
    grid_BOLD(:,:,i) = F(mesh_xcoords, mesh_ycoords);
end

grid_BOLD(isnan(grid_BOLD)) = 0;

%% Saving grid_BOLD and F in a mat file

if saveMATfile
    filename = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/',hemisphere,...
                '.Scan',num2str(scanNo),'_resolution=',num2str(resolution),'.mat'];
    save(filename, 'grid_BOLD', 'F')
end
