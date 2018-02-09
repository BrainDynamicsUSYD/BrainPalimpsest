function [templates, eccentricity, polar, responses_1D] = ...
    calculatingV1AverageTimeEvolution(hemisphere, resolution, ...
    polar_deviation, avg_BOLD_signal, deconvResponses_avg, isGridBenson, saveGridBenson)
%% calculatingV1AverageTimeEvolution.m
%
% Calculates the average of the 1D responses in V1 and converts them in 
% 1D form (eccentricity, time).
%
% Inputs: hemisphere        : string of hemisphere
%                             Possible fields are lh for left hemisphere and 
%                             rh for right hemisphere.
%         resolution        : x and y spatial resolution in mm
%         polar_deviation   : polar angle deviation from 90deg that will be 
%                             considered for the averaging process
%         avg_BOLD_signal   : array of averaged BOLD signal (x,y,t)
%                             size(avg_BOLD_signal) = [length(y), length(x), length(t)]
%         deconvResponses_avg    : array of averaged deconvolved responses (x,y,t)
%                             size(deconvResponses_avg.{}) = [length(y), length(x), length(t)]
%                             Possible fields are reconvBOLD, neural, 
%                             neuroglial, CBF, CBV, dHb, Wmode, Lmode, and Dmode.
%         isGridBenson      : 1 or 0
%                             Choose 1 if a gridded form of the Benson 
%                             map is available, 0 otherwise.
%         saveGridBenson    : 1 or 0
%                            Choose 1 if you want to save the mat file of 
%                            the gridded Benson map, 0 otherwise.
%
% Output: templates         : struct of area, polar_angle, eccentricity 
%                             templates from the Benson map
%                             Possible fields are area, polar_angle, eccentricity, 
%                             area_v1, polar_v1, and eccentricity_v1.
%         eccentricity      : struct
%                             Possible fields are values, rows, cols, and centroids.  
%         polar             : struct
%                             Possible fields are deviations, rows, and cols.  
%         responses_1D      : array of processed 1D responses in V1
%                             Possible fields are BOLD, reconvBOLD, neural, neuroglial, 
%                             CBF, CBV, dHb, Wmode, Lmode, and Dmode.
% 
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%% Loading relevant files

if isGridBenson
    filename_Benson = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/',hemisphere,...
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

%% Finding voxels within the same eccentricity bin

eccentricity_diff = 0.05;
eccentricity.values = eccentricity_diff:eccentricity_diff:5+eccentricity_diff;

for i = 1:length(eccentricity.values)-1
    [row, col] = find((templates.eccentricity_v1 >= eccentricity.values(i) & ...
                   templates.eccentricity_v1 < eccentricity.values(i+1)));
    if strcmp(hemisphere, 'lh')
        ind = find(row > 210);
    elseif strcmp(hemisphere, 'rh')
        ind = find(row < 270);
    end
    row = row(ind);  
    col = col(ind);
    eccentricity.rows{i} = row';
    eccentricity.cols{i} = col';
    eccentricity.centroids(i, :) = [round(mean(col)), round(mean(row))];
end

%% Finding voxels within the same polar angle bin

polar.deviations = 1:81;

for i = 1:length(polar.deviations)-1
    [row, col] = find((templates.polar_v1 > 90-polar.deviations(i) & templates.polar_v1 < 90) | ...
                      (templates.polar_v1 > 90 & templates.polar_v1 < 90+polar.deviations(i)));
    if strcmp(hemisphere, 'lh')
        ind = find(row > 210);
    elseif strcmp(hemisphere, 'rh')
        ind = find(row < 270);
    end
    row = row(ind);
    col = col(ind);
    polar.rows{i} = row';
    polar.cols{i} = col';
end

%% Calculating 1D responses as a function of time and eccentricity

responses = {'BOLD', 'reconvBOLD', 'neural', 'neuroglial', 'CBF', 'CBV', ...
             'dHb', 'Wmode', 'Lmode', 'Dmode'};

for i = 1:length(responses)
    if i==1
        data = real(avg_BOLD_signal);
    else
        data = real(deconvResponses_avg.(responses{i}));
    end
    
    data(isnan(data)) = 0;
    
    data_1D = zeros(size(data, 3), length(eccentricity.values)-1);
    
    for j = 1:size(data_1D, 1)
        for k = 1:size(data_1D, 2)
            ecc = [eccentricity.rows{k}', eccentricity.cols{k}'];
            pol = [polar.rows{polar_deviation}', polar.cols{polar_deviation}'];
            intersection = intersect(ecc, pol, 'rows');
            
            data_1D(j, k) = mean2(data(intersection(:,1), ...
                               intersection(:,2), j));
        end
    end
    
    responses_1D.(responses{i}) = data_1D;
end
