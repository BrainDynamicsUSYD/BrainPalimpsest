%% Results_ExpandingRingAndExpandingArc.m
%
% This script shows all the necessary steps to produce the results for 
% both the 2D Expanding Ring and Expanding Arc data of the manuscript. 
%
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%% Adding the paths of the sub-directories for direct access of files
% This is not necessary if entire Palimpsest toolbox is added to the Matlab 
% path via addpath(genpath('PalimpsestToolboxLocation')) where
% PalimpsestToolboxLocation is the location of the toolbox

addpath('Data', 'Functions', 'PlottingFunctions')

%% Loading the default values of the model parameters 

params = loadParameters;

%% Loading the data

UseSampleData = 1;

if UseSampleData   
    % Load sample experimental expanding ring or expanding arc
    
    % Set which hemisphere to analyze
    % lh for left hemisphere and rh for right hemisphere
    hemisphere = 'lh'; 
    
    % Set which data to be loaded
    % either 'ExpandingRing' or 'ExpandingArc'
    DataToDemonstrate = 'ExpandingRing';
%     DataToDemonstrate = 'ExpandingArc';
    
    % Set which scan number of the data to be loaded
    % scanNo = 1 for the expanding ring data
    % scanNo = 2, 3, ..., 11 for the expanding arc data
    if strcmpi(DataToDemonstrate, 'ExpandingRing')
        scanNo = 1;
    elseif strcmpi(DataToDemonstrate, 'ExpandingArc')
        scanNo = 2;             % change according to which scan number 
                                % (2 to 11) you want to see
    end
    
    % If flattened gridded (matrix) form of data is available, set isGridBOLD = 1
    % If data is in flattened triangulated form and gridded form is not available, 
    % set isGridBOLD = 0 and the data will be converted to gridded form
    isGridBOLD = 0;
    
    % Set x and y spatial resolution in mm of the gridded data
    resolution = 0.2;
    
    % load the actual data
    % change according to your data
    % 3D matrix must be oriented in this way: 
    % rows = y, columns = x, third = time
    if isGridBOLD
        filename = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/',hemisphere,...
                    '.Scan',num2str(scanNo),'_resolution=',num2str(resolution),'.mat'];
        load(filename, 'grid_BOLD', 'F')        % grid_BOLD is the data in 
                                                % matrix form and F is the 
                                                % interpolant object 
    else
        [grid_BOLD, F] = makingGriddedFlat_BOLD(hemisphere, scanNo, ...
                                                resolution, 0);
    end
    
    BOLD_signal = grid_BOLD;
    avg_BOLD_signal = averagingTimeSeries(scanNo, BOLD_signal);
else
    % If have your own data, load it here
    % Follow the above format and peform some necessary changes to the
    % relevant functions.
end

%% Loading experimental times and distances

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

xcoords = flat_struct.Vertices(:, 1);
ycoords = flat_struct.Vertices(:, 2);

x_lim = floor(max(xcoords));
y_lim = floor(max(ycoords));

% experimental distances and times
x_experiment = (-x_lim:resolution:x_lim)*1e-3; % need to multiply by 1e-3 to convert to m
y_experiment = (-y_lim:resolution:y_lim)*1e-3; % need to multiply by 1e-3 to convert to m

dt = 2;
t_experiment_orig = (0:size(BOLD_signal, 3)-1)*dt;      % dt interval
t_experiment_avg = (0:size(avg_BOLD_signal, 3)-1)*dt;   % dt interval

%% Choosing which astrocytic model scheme will be used
% 'orig_noast' : original model values for wf and kappa + without astrocytic delay
% 'orig_ast'   : original model values for wf and kappa + with astrocytic delay
% 'new_noast'  : new model values for wf and kappa + without astrocytic delay
% 'new_ast'    : new model values for wf and kappa + with astrocytic delay
% 'mean_ast'   : mean experimental values for wf, kappa, and astrocytic delay

astrocyte_scheme = 'mean_ast';

% Mean estimates of wf, kappa, and astrocytic delay will be used depending on
% which astrocytic model scheme is used
load Data/ExpandingRingAndExpandingArc/parameterEstimates_AstrocyticDelay.mat

% set wf, kappa, and tau_d to the experimental mean values 
params.w_f = results.(astrocyte_scheme).w_f;
params.kappa = results.(astrocyte_scheme).kappa;
params.tau_d = results.(astrocyte_scheme).tau_d;

%% Setting some processing info
polar_deviation = 20;       % polar angle deviation from 90deg that will be 
                            % considered for the averaging of the time 
                            % evolution of responses in V1
isGridBenson = 0;           % 1 if a gridded form of the Benson map is 
                            % available, 0 otherwise
isVisualStimulus = 0;       % 1 if visual stimuli were already created,
                            % 0 otherwise

saveGridBenson = 1;         % save gridded Benson map in a mat file
saveMAT = 1;                % save deconvolution results in a mat file
saveOverlay = 1;            % save deconvolution results in freesurfer compatible file
                            % overlay format in a mgz file

%% Processing the BOLD signal by refining resolution and padding zeros

disp('Processing resolution of BOLD signal ...')

params.Nkx = 2^7;
params.Nky = 2^7;
params.Nw = 2^8;

% BOLD_processed is the processed BOLD signal
% x, y, and t are the new distance and time vectors, respectively
[BOLD_processed, x, y, t] = BOLD_processing_2D(BOLD_signal, x_experiment, ...
                                       y_experiment, t_experiment_orig, params);
                                
%% Deconvolution of 2D responses 

disp('Deconvolving the responses ...')

% changing the NSR term
params.noise = 0.5;             % constant

deconvResponses = deconvolution_HybridWiener_2D(BOLD_processed, x, y, t, params);

%% Return responses back to experimental x, y, and t space

disp('Processing the responses to return to experimental x, y, t ...')

responses = {'BOLD', 'reconvBOLD', 'neural', 'neuroglial', 'CBF', 'CBV', ...
             'dHb', 'Wmode', 'Lmode', 'Dmode'};

x_index = dsearchn(x', x_experiment')';
y_index = dsearchn(y', y_experiment')';
t_index = dsearchn(t', t_experiment_orig')';
    
for k = 2:length(responses)
    data = real(deconvResponses.(responses{k}));
    
    deconvResponses.(responses{k}) = data(y_index, x_index, t_index);
    deconvResponses_avg.(responses{k}) = averagingTimeSeries(scanNo, data(y_index, x_index, t_index));
end

%% Calculating the average time evolution of different eccentricities in V1

disp('Calculating 1D responses in V1...')

[templates, eccentricity, polar_angle, responses_1D] = ...
    calculatingV1AverageTimeEvolution(hemisphere, resolution, polar_deviation, ...
    avg_BOLD_signal, deconvResponses_avg, isGridBenson, saveGridBenson);

%% Reordering the time dimension of gridded avgBOLD, gridded avg deconvResponses, 
%  and avg 1D responses
% This block is only relevant for the expanding ring data

if strcmpi(DataToDemonstrate, 'ExpandingRing')
    t_reorder_index = [8:15, 1:7];
    disp('Reordering the time dimension of average results ...')

    for k = 1:length(responses)
        if k == 1
            reordered_avg_BOLD_signal = avg_BOLD_signal(:, :, t_reorder_index);
        else
            data1 = real(deconvResponses_avg.(responses{k}));
            reordered_deconvResponses_avg.(responses{k}) = data1(:, :, t_reorder_index);
        end
        data2 = real(responses_1D.(responses{k}));
        reordered_responses_1D.(responses{k}) = data2(t_reorder_index, :);
    end
end

%% Making the visual stimuli of the data in gridded form 

if UseSampleData
    disp('Making visual stimuli ...')
    
    if isVisualStimulus
        filename = ['Data/ExpandingRingAndExpandingArc/VisualStimulus/', ...
                    hemisphere,'.Scan',num2str(scanNo),...
                    '_VisualStimulus_resolution=',num2str(resolution),'.mat'];
        load(filename, 'v1_boundary', 'thmat', 'thmat_templateSpace', 'rmat', ...
                   'visualStimulus_raw', 'visualStimulus_smooth')
    else
        num_stimuli = length(t_experiment_avg);

        [v1_boundary, thmat, thmat_templateSpace, rmat, visualStimulus_raw, ...
            visualStimulus_smooth] = makingVisualStimulus(hemisphere, resolution, ...
                        DataToDemonstrate, num_stimuli, isGridBenson, saveGridBenson);
    end
end
 
%% Saving the gridded BOLD, gridded avgBOLD, gridded deconvResponses, 
%   gridded avg deconvResponses, 1D responses, reordered avgBOLD, 
%   reordered avg deconvResponses, reordered 1D responses, x, y, t, F, Benson 
%   templates, eccentricity and polar angle values of 1D responses, 
%   v1_boundary, thmat, thmat_templateSpace, rmat, visualStimulus_raw, and
%   visualStimulus_smooth on a mat file for future use

if saveMAT
    disp('Saving results in a mat file ...')
    
    filename = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/',hemisphere,...
                '.Scan',num2str(scanNo),'_resolution=',num2str(resolution),'.mat'];
            
    if strcmpi(DataToDemonstrate, 'ExpandingRing')
        save(filename, 'grid_BOLD', 'avg_BOLD_signal', 'deconvResponses', ...
                       'deconvResponses_avg', 'responses_1D', 'reordered_avg_BOLD_signal', ...
                       'reordered_deconvResponses_avg', 'reordered_responses_1D', ...
                       'x_experiment', 'y_experiment', 't_experiment_orig', ...
                       't_experiment_avg', 'F', 'templates', 'eccentricity', ... 
                       'polar_angle')
    elseif strcmpi(DataToDemonstrate, 'ExpandingArc')
        save(filename, 'grid_BOLD', 'avg_BOLD_signal', 'deconvResponses', ...
                       'deconvResponses_avg', 'responses_1D', ...
                       'x_experiment', 'y_experiment', 't_experiment_orig', ...
                       't_experiment_avg', 'F', 'templates', 'eccentricity', ...
                       'polar_angle')
    end
    
    if ~isVisualStimulus
        filename = ['Data/ExpandingRingAndExpandingArc/VisualStimulus/', ...
                        hemisphere,'.Scan',num2str(scanNo),...
                        '_VisualStimulus_resolution=',num2str(resolution),'.mat'];
        save(filename, 'v1_boundary', 'thmat', 'thmat_templateSpace', 'rmat', ...
                       'visualStimulus_raw', 'visualStimulus_smooth')
    end
end

%% Saving the overlay files for freesurfer

if saveOverlay
    disp('Saving overlay results in a mgz file ...')
    
    for i = 1:length(responses)
        if i==1
            data_orig = BOLD_signal;
            data_avg = avg_BOLD_signal;
            if strcmpi(DataToDemonstrate, 'ExpandingRing')
                reordered_data_avg = reordered_avg_BOLD_signal;
            end
        else
            data_orig = deconvResponses.(responses{i});
            data_avg = deconvResponses_avg.(responses{i});
            if strcmpi(DataToDemonstrate, 'ExpandingRing')
                reordered_data_avg = reordered_deconvResponses_avg.(responses{i});
            end
        end

        overlay_orig = goingBackToTriangulation(hemisphere, scanNo, ...
                         data_orig, F, x_experiment, y_experiment);
        MRIwrite(overlay_orig, ['Data/ExpandingRingAndExpandingArc/ResponseSurfaces/', ...
            hemisphere,'.Scan',num2str(scanNo),'_',responses{i},'_orig.mgz']);

        overlay_avg = goingBackToTriangulation(hemisphere, scanNo, ...
                         data_avg, F, x_experiment, y_experiment);
        MRIwrite(overlay_avg, ['Data/ExpandingRingAndExpandingArc/ResponseSurfaces/', ...
            hemisphere,'.Scan',num2str(scanNo),'_',responses{i},'_avg.mgz']);
        
        if strcmpi(DataToDemonstrate, 'ExpandingRing')
            overlay_reordered_avg = goingBackToTriangulation(hemisphere, scanNo, ...
                             reordered_data_avg, F, x_experiment, y_experiment);
            MRIwrite(overlay_reordered_avg, ['Data/ExpandingRingAndExpandingArc/ResponseSurfaces/', ...
                hemisphere,'.Scan',num2str(scanNo),'_',responses{i},'_reordered_avg.mgz']);
        end
    end
    
    % For visual stimulus
    overlay_visualStimulus_raw = goingBackToTriangulation(hemisphere, scanNo, ...
                         visualStimulus_raw, F, x_experiment, y_experiment);
    MRIwrite(overlay_visualStimulus_raw, ['Data/ExpandingRingAndExpandingArc/VisualStimulus/', ...
             hemisphere,'.Scan',num2str(scanNo),'_VisualStimulus_raw.mgz']);

    overlay_visualStimulus_smooth = goingBackToTriangulation(hemisphere, scanNo, ...
                     visualStimulus_smooth, F, x_experiment, y_experiment);
    MRIwrite(overlay_visualStimulus_smooth, ['Data/ExpandingRingAndExpandingArc/VisualStimulus/', ...
             hemisphere,'.Scan',num2str(scanNo),'_VisualStimulus_smooth.mgz']);
end

%% Plotting the result 

disp('Plotting and saving the results ...')

% boolean variable
% choose normalization = 1 if you want the responses to be normalized, 0
% otherwise
% Note that normalization is different for each quantity.
%   1. Since we do not have a baseline neural activity, we normalize the
%   response with respect to maximum.
%   2. Since we have baseline CBF, CBV, and dHb, we normalize their
%   responses with respect to the baseline values.
%   3. Since the modes do not have baseline values, we normalize their
%   responses with respect to the maximum of the total BOLD response.
normalization = 0;

% zoom_lim changes the limits of the patch of the cortex to be shown 
% code below will show the entire coverage of the experimental data
zoom_lim = [min(x_experiment)*1e3, max(x_experiment)*1e3;
            min(y_experiment)*1e3, max(y_experiment)*1e3];

% string pertaining to which response will be plotted
% Possible CASE-INSENSITIVE inputs are: 
% 'BOLD', 'reconvBOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', 
% 'Wmode', 'Lmode', 'Dmode', 'all_w_BOLD', 'all_no_BOLD' 

% Note that 'all_w_BOLD' combines all the responses (except reconvBOLD) + BOLD 
% in a single figure, 
% while 'all_no_BOLD' combines all the responses (except reconvBOLD) without 
% BOLD in a single figure
plot_what = 'all_w_BOLD';

% string pertaining to which quantity will the other quantities be
% correlated with
% Possible CASE-INSENSITIVE inputs are: 
% 'visualStimulus', 'BOLD', 'reconvBOLD', 'neural', 'neuroglial', 'CBF', 
% 'CBV', 'dHb', 'Wmode', 'Lmode', 'Dmode'
what_correlation = 'BOLD';


if strcmpi(DataToDemonstrate, 'ExpandingRing')
    % Time slices of responses in (x, y, t = tslice)
    tslice = 0:4:16;                % time slices of responses to show
    clim_factor1 = 5;               % decreasing the limits of the colorbar
    fig1 = plot2D_VisualMap(visualStimulus_smooth, x_experiment, y_experiment, t_experiment_avg, ...
                    tslice, zoom_lim, v1_boundary);
    fig2 = plot2D_TimeSlices(BOLD_signal, deconvResponses, x_experiment, y_experiment, t_experiment_orig, ...
                    tslice, zoom_lim, params, normalization, clim_factor1, plot_what, v1_boundary);
    fig3 = plot2D_TimeSlices(reordered_avg_BOLD_signal, reordered_deconvResponses_avg, ...
                             x_experiment, y_experiment, t_experiment_avg, ...
                             tslice, zoom_lim, params, normalization, clim_factor1, plot_what, v1_boundary);
                

    % Evolution of responses as a movie                        
    movie_filename = ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '.avi'];
    frame_rate = mean(diff(t))*2;                   % frames per second of movie
    
    fig4 = plot2D_Movie(reordered_avg_BOLD_signal, reordered_deconvResponses_avg, x_experiment, y_experiment, t_experiment_avg, ...
                 zoom_lim, params, normalization, clim_factor1, plot_what, movie_filename, frame_rate, v1_boundary);
    
    
    % Time series of V1 responses in (eccentricity, t)
    clim_factor2 = 1;
    fig5 = plot_V1TimeEvolution(reordered_responses_1D, eccentricity.values(1:end-1), ...
                                t_experiment_avg, params, normalization, clim_factor2, plot_what);
                            
    
    % Normalized time profiles of V1 responses for different eccentricities
    eccentricity_interest = [0.5, 1.5, 2.5];
    fig6 = plot_V1TimeProfilePerEccentricity(reordered_responses_1D, eccentricity.values(1:end-1), ...
                                        eccentricity_interest, t_experiment_avg, plot_what);
                                    
    % Normalized time profiles of V1 responses for one eccentricity
    eccentricity_interest = 1.5;
    fig7 = plot_V1TimeProfilePerEccentricity(reordered_responses_1D, eccentricity.values(1:end-1), ...
                                        eccentricity_interest, t_experiment_avg, plot_what);
    
elseif strcmpi(DataToDemonstrate, 'ExpandingArc')
    % Time slices of responses in (x, y, t = tslice)
    tslice = 0:4:16;                % time slices of responses to show
    clim_factor1 = 5;               % decreasing the limits of the colorbar
    fig1 = plot2D_VisualMap(visualStimulus_smooth, x_experiment, y_experiment, t_experiment_avg, ...
                    tslice, zoom_lim, v1_boundary);
    fig2 = plot2D_TimeSlices(BOLD_signal, deconvResponses, x_experiment, y_experiment, t_experiment_orig, ...
                    tslice, zoom_lim, params, normalization, clim_factor1, plot_what, v1_boundary);
    fig3 = plot2D_TimeSlices(avg_BOLD_signal, deconvResponses_avg, ...
                             x_experiment, y_experiment, t_experiment_avg, ...
                             tslice, zoom_lim, params, normalization, clim_factor1, plot_what, v1_boundary);

                            
    % Evolution of responses as a movie                        
    movie_filename = ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '.avi'];
    frame_rate = mean(diff(t))*2;                   % frames per second of movie
    
    fig4 = plot2D_Movie(avg_BOLD_signal, deconvResponses_avg, x_experiment, y_experiment, t_experiment_avg, ...
                 zoom_lim, params, normalization, clim_factor1, plot_what, movie_filename, frame_rate, v1_boundary);
    
    % Time series of V1 responses in (eccentricity, t)
    clim_factor2 = 1;
    fig5 = plot_V1TimeEvolution(responses_1D, eccentricity.values(1:end-1), ...
                                t_experiment_avg, params, normalization, clim_factor2, plot_what);
                            
                            
    % Normalized time profiles of V1 responses for different eccentricities
    eccentricity_interest = [0.5, 1.5, 2.5];
    fig6 = plot_V1TimeProfilePerEccentricity(responses_1D, eccentricity.values(1:end-1), ...
                                        eccentricity_interest, t_experiment_avg, plot_what);
                                    
    % Normalized time profiles of V1 responses for one eccentricity
    eccentricity_interest = 1.5;
    fig7 = plot_V1TimeProfilePerEccentricity(responses_1D, eccentricity.values(1:end-1), ...
                                        eccentricity_interest, t_experiment_avg, plot_what);
end

% Cross-correlations of the response corresponding to what_correlation
% with other responses
fig8 = plot_CrossCorrelations(hemisphere, what_correlation);
    
set(fig1, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig2, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig3, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig5, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig6, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig7, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig8, 'PaperPositionMode','auto')     %# WYSIWYG
print(fig1, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_VisualStimulus.eps'])
print(fig2, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_TimeSlices.eps'])
print(fig3, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_TimeSlices_avg.eps'])
print(fig5, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_V1TimeEvolution.eps'])
print(fig6, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_V1TimeProfilePerEccentricity.eps'])
print(fig7, '-painters', '-depsc', ['Figures/ExperimentalResults_', DataToDemonstrate, '_Scan',num2str(scanNo), '_V1TimeProfileOneEccentricity.eps'])
print(fig8, '-painters', '-depsc', ['Figures/ExperimentalResults_ExpandingRingAndExpandingArc_CrossCorrelationWith_', what_correlation, '.eps'])

%% END

close(fig1); close(fig2); close(fig3); 
close(fig4); close(fig5); close(fig6);
close(fig7); close(fig8);
disp('Finished')
disp([])
