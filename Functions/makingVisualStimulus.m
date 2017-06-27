function [v1_boundary, thmat, thmat_templateSpace, rmat, visualStimulus_raw, ...
    visualStimulus_smooth] = makingVisualStimulus(hemisphere, resolution, ...
    DataToDemonstrate, num_stimuli, isGridBenson, saveGridBenson)
%% makingVisualStimulus.m
%
% Makes the visual stimulus of the expanding ring and expanding wedge data
% in grid form using the Benson template.
%
% Inputs: hemisphere        : string of hemisphere
%                             Possible fields are lh for left hemisphere and 
%                             rh for right hemisphere.
%         resolution        : x and y spatial resolution in mm
%         DataToDemonstrate : string
%                             Possible fields are 'ExpandingRing' and 
%                             'ExpandingWedge'.
%         num_stimuli       : number of stimuli to be produced 
%         isGridBenson      : 1 or 0
%                             Choose 1 if a gridded form of the Benson 
%                             map is available, 0 otherwise.
%         saveGridBenson    : 1 or 0
%                            Choose 1 if you want to save the mat file of 
%                            the gridded Benson map, 0 otherwise.
%
% Output: v1_boundary       : structure of boundary of v1
%                             Possible fields are x for the x pixel
%                             coordinate, y for the y pixel coordinate, and
%                             matrix for its matrix form
%         thmat             : matrix of polar angle mappings of the visual stimulus
%         thmat             : matrix of polar angle mappings of the visual stimulus
%                             in Benson template space
%         rmat              : matrix of eccentricity mappings of the visual stimulus
%         visualStimulus_raw    : matrix of actual visual stimulus
%         visualStimulus_smooth : matrix of smoothed version of the visual
%                                 stimulus
% 
% James Pang, University of Sydney, 2017

%% Defining necessary parameters

a                   =   0.75;             % Foveal pole
K                   =   18;               % scaling parameter
maxScreenEc         =   5.5;              % Maxium eccentricity of the screen in the y-direction.
refinement_stim     =   500;              % refinement of the stimulus, the higher this is the more points there will be
bar_width           =   1;                % width of the bar (in mm)
bar_height          =   10;               % height of the bar (in mm)
maxCortex           =   32;               % maximum distance on canonical cortex the stimulus will travel
d_s                 = (maxCortex-1)/(num_stimuli-1); % stepsize of the stimulus, the end point will reach maxCortex size

%% Loading relevant files

if isGridBenson
    filename_Benson = ['Data/ExpandingRingAndExpandingWedge/GriddedMatFiles/',hemisphere,...
                       '.Benson_resolution=',num2str(resolution),'.mat'];
    load(filename_Benson, 'grid_Benson')
else
    [grid_Benson, ~] = makingGriddedFlat_Benson(hemisphere, resolution, ...
                                            saveGridBenson);
end
              
%% Defining the area, polar, and eccentricity templates

templates.polar_angle = grid_Benson(:,:,1);
templates.eccentricity = grid_Benson(:,:,2);
templates.area = grid_Benson(:,:,3);

templates.area_v1 = (templates.area >= -1.0001 & templates.area <= -0.0001) ...
                    | (templates.area <= 1.0001 & templates.area >= 0.0001);
templates.polar_v1 = templates.polar_angle.*templates.area_v1;
templates.eccentricity_v1 = templates.eccentricity.*templates.area_v1;

%% Defining the boundary of v1

v1 = templates.area_v1;
if strcmp(hemisphere, 'lh')
    v1(1:210,:) = 0;
elseif strcmp(hemisphere, 'rh')
    v1(270:end,:) = 0;
end
B = bwboundaries(v1, 4, 'noholes');
v1_boundary.x = B{1}(:,2);
v1_boundary.y = B{1}(:,1);
v1_boundary.matrix = zeros(size(v1));
for i=1:length(v1_boundary.x);
    v1_boundary.matrix(v1_boundary.y(i), v1_boundary.x(i)) = 1;
end

%% Making bar stimulus

x = zeros(num_stimuli, refinement_stim);
y = zeros(num_stimuli, refinement_stim);
for k=1:num_stimuli
    x(k,:) = linspace(1, bar_width, refinement_stim) + (k-1)*d_s;
    y(k,:) = linspace(-bar_height/2, bar_height/2, refinement_stim); 
end

%% Converting to polar coordinates (eccentricity and polar angle)

for k=1:num_stimuli
    [xx, yy] = meshgrid(x(k,:), y(k,:));  
    [TH, R] = cart2pol(xx, yy);
    
    w = exp(1i*reshape(TH, numel(TH), 1)).*reshape(R, numel(R), 1);
    z = exp(w/K) - a;
    
    if strcmpi(DataToDemonstrate, 'ExpandingRing')
        thVF = linspace(-pi/2, pi/2, length(x(k,:)))';
        rVF = mean(abs(z))*ones(size(thVF));
        combined = [thVF, rVF];
        unique_vals = unique(combined, 'rows', 'stable');
    elseif strcmpi(DataToDemonstrate, 'ExpandingWedge')
        thVF = angle(z);
        rVF = abs(z);
        combined = [thVF, rVF];
        unique_vals = unique(combined, 'rows', 'stable');
    end
    
    if strcmp(hemisphere, 'lh')
        thmat(:,k) = unique_vals(:,1);
        thmat_templateSpace(:, k) = rad2deg(unique_vals(:,1) + pi/2); 
        rmat(:,k) = unique_vals(:,2);
    elseif strcmp(hemisphere, 'rh')
        thmat(:,k) = unique_vals(:,1) + pi;
        thmat_templateSpace(:, k) = rad2deg(unique_vals(:,1) + pi - pi/2); 
        rmat(:,k) = unique_vals(:,2);
    end
end

%% Projecting in a flattened Benson template
        
ecc_thresh = 0.1;       % eccentricity threshold 
polar_thresh = 1;       % polar angle threshold

visualStimulus_raw = zeros(size(templates.eccentricity, 1), ...
                           size(templates.eccentricity, 2), ...
                           num_stimuli);
visualStimulus_smooth = zeros(size(templates.eccentricity, 1), ...
                           size(templates.eccentricity, 2), ...
                           num_stimuli);
for k=1:num_stimuli
    intersections = [];
    for j=1:size(rmat,1)
        [ecc_row, ecc_col] = find(abs(templates.eccentricity_v1 - rmat(j,k)) <= ecc_thresh);
        [polar_row, polar_col] = find(abs(templates.polar_v1 - thmat_templateSpace(j,k)) <= polar_thresh);
        
        if strcmp(hemisphere, 'lh')
            ecc_ind_v1 = find(ecc_row > 210);
            polar_ind_v1 = find(polar_row > 210);
        elseif strcmp(hemisphere, 'rh')
            ecc_ind_v1 = find(ecc_row < 270);
            polar_ind_v1 = find(polar_row < 270);
        end
        
        ecc_locations = [ecc_row(ecc_ind_v1), ecc_col(ecc_ind_v1)];
        polar_locations = [polar_row(polar_ind_v1), polar_col(polar_ind_v1)];

        intersections = [intersections; ecc_locations(ismember(ecc_locations, ...
                                              polar_locations, 'rows'), :)];
        
        if ~isempty(intersections)
            visualStimulus_raw(intersections(end,1), intersections(end,2),k) = 1;
        end
    end
    
    kernel = ones(5);    % smoothing kernel
    smooth = conv2(visualStimulus_raw(:,:,k), kernel, 'same'); 
    smooth(smooth~=0) = 1;
    visualStimulus_smooth(:,:,k) = smooth;
end

% Visualization
%
% figure;
% for k=1:num_stimuli
%     p2 = polar(NaN, 5.5, '.');
%     hold on;
%     p3 = polar(thmat(:,k), rmat(:,k), '.');
%     hold off;
%     pause(0.5)
% end
% 
% figure;
% for k=1:num_stimuli
%     subplot(1,2,1)
%     imagesc(visualStimulus_raw(:,:,k))
%     set(gca,'ydir','normal')
%     colormap(flipud(bone))
%     axis image
%     
%     subplot(1,2,2)
%     imagesc(visualStimulus_smooth(:,:,k))
%     set(gca,'ydir','normal')
%     colormap(flipud(bone))
%     axis image
%     
%     pause(0.5)
% end