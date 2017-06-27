%% Results_StationaryRing.m
%
% This script shows all the necessary steps to produce the results for 
% the 1D Stationary Ring data of the manuscript. 
%
% James Pang, University of Sydney, 2016

%% Adding the paths of the sub-directories for direct access of files
% This is not necessary if entire Palimpsest toolbox is added to the Matlab 
% path via addpath(genpath('PalimpsestToolboxLocation')) where
% PalimpsestToolboxLocation is the location of the toolbox

% addpath('Data', 'Functions', 'PlottingFunctions')

%% Loading the default values of the model parameters 

params = loadParameters;
 
%% Loading the data

UseSampleData = 1;

if UseSampleData   
    % Load sample experimental data from Aquino et al. (PLoS Comp Biol, 2012)
    
    % Set subject name 
    % Possible subject names are 'aLeft', 'aRight', 'rLeft', 'rRight', 
    % 'sLeft', 'sRight', 'vRight'.
    subject = 'aRight';
    
    % Get the estimated v_b and Gamma for subject
    load Data/StationaryRing/fileParams.mat

    % Find index of subject in associatedParams cell
    index = find(strcmp(associatedParams, subject));
    
    % Assigns the estimated experimental v_b and Gamma to variables for 
    % future use
    params.v_b = Vmat(index)*1e-3;      % Vmat is in mm/s so multiply by 1e-3 to convert to m/s
    params.Gamma = -Gmat(index);        % Gmat is written as -Gamma so extra - sign is need
    
    % load the actual data
    % change according to your data
    % matrix must be oriented in this way: 
    % rows = x, columns = time
    load(cat(2, 'Data/StationaryRing/', subject, '.mat'))

    BOLD_signal = smoothTSER;           % actual 1D BOLD data from toolbox

    % experimental times and distances from Aquino et al. (PLoS Comp Biol, 2012)
    % change according to your protocol
    t_experiment = 0.25:0.25:20;        % 20 seconds with 250ms interval
    x_experiment = distanceInterp*1e-3; % need to multiply by 1e-3 to convert to m
else
    % If have your own data, load it here
    % v_b and Gamma will be set to their default values of 2 mm/s and 0.8/s 
    % based on experimental mean measurements
    % Note that you shall follow the above steps such as providing
    % t_experiment and x_experiment according to your protocol
end

%% Choosing which astrocytic model scheme will be used
% 'orig_noast' : original model values for wf and kappa + without astrocytic delay
% 'orig_ast'   : original model values for wf and kappa + with astrocytic delay
% 'new_noast'  : new model values for wf and kappa + without astrocytic delay
% 'new_ast'    : new model values for wf and kappa + with astrocytic delay
% 'mean_ast'   : mean experimental values for wf, kappa, and astrocytic delay

astrocyte_scheme = 'new_ast';

if UseSampleData
    % load astrocytic parameter estimates. 
    load Data/StationaryRing/parameterEstimates_AstrocyticDelay.mat
    
    % set estimated values of wf, kappa, and astrocytic delay tau_d
    params.w_f = results.(astrocyte_scheme).w_f(index);
    params.kappa = results.(astrocyte_scheme).kappa(index);
    params.tau_d = results.(astrocyte_scheme).tau_d(index);
else
    % If you have your own data, mean estimates of wf, kappa, and astrocytic 
    % delay will be used 
    
    astrocyte_scheme = 'mean_ast';
    
    % set wf, kappa, and tau_d to the experimental mean values 
    params.w_f = results.(astrocyte_scheme).w_f;
    params.kappa = results.(astrocyte_scheme).kappa;
    params.tau_d = results.(astrocyte_scheme).tau_d;
end

%% Processing the BOLD signal by refining resolution and padding zeros

display('Processing resolution of BOLD signal ...')

% BOLD_processed is the processed BOLD signal
% x and t are the new distance and time vectors, respectively
[BOLD_processed, x, t] = BOLD_processing_1D(BOLD_signal, ...
                                    x_experiment, t_experiment, params);
                                
%% Deconvolution of 1D responses 

display('Deconvolving the responses ...')

% changing the NSR term
params.noise = 0.5;             % constant
% params.noise = [500, 0.1];    % cut-off frequencies

deconvResponses = deconvolution_HybridWiener_1D(BOLD_processed, x, t, params);

%% Plotting the result 

display('Plotting and saving the results ...')

% boolean variable
% Choose normalization = 1 if you want the responses to be normalized, 0
% otherwise
% Note that normalization is different for each quantity.
%   1. Since we do not have a baseline neural activity, we normalize the
%   response with respect to maximum.
%   2. Since we have baseline CBF, CBV, and dHb, we normalize their
%   responses with respect to their baseline values.
%   3. Since the modes do not have baseline values, we normalize their
%   responses with respect to the maximum of the total BOLD signal.
normalization = 0; 


% string pertaining to which response(s) will be plotted
% Possible CASE-INSENSITIVE inputs are: 
% 'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 'Dmode', 
% 'all_w_BOLD', 'all_no_BOLD' 

% Note that 'all_w_BOLD' combines all the responses + BOLD in a single figure,
% while 'all_no_BOLD' combines all the responses without BOLD in a single figure
plot_what = 'all_w_BOLD';


% Contour maps of responses in (x, t)
fig1 = plot1D_SpatioTemporal(BOLD_processed, deconvResponses, x, t, params, ...
                             normalization, plot_what);

% Time series of responses in (x=0, t)
fig2 = plot1D_CenterTimeSeries(BOLD_processed, deconvResponses, x, t, params, ...
                             normalization, plot_what);

% Combined contour maps and time series at x=0 of responses
fig3 = plot1D_Combined_SpatioTemporal_CenterTimeSeries(BOLD_processed, ...
                deconvResponses, x, t, params, normalization, plot_what);

% Evolution of responses as a movie
movie_filename = 'Figures/ExperimentalResults_StationaryRing.avi';
frame_rate = 10;        % frames per second of movie

fig4 = plot1D_Movie(BOLD_processed, deconvResponses, x, t, params, ...
             normalization, plot_what, movie_filename, frame_rate);
         
set(fig1, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig2, 'PaperPositionMode','auto')     %# WYSIWYG
set(fig3, 'PaperPositionMode','auto')     %# WYSIWYG
print(fig1, '-painters', '-depsc', 'Figures/ExperimentalResults_StationaryRing_SpatioTemporal.eps')
print(fig2, '-painters', '-depsc', 'Figures/ExperimentalResults_StationaryRing_CenterTimeSeries.eps')
print(fig3, '-painters', '-depsc', 'Figures/ExperimentalResults_StationaryRing_CombinedSpatiotemporalAndCenterTimeSeries.eps')
      
%% END

close(fig1); close(fig2); close(fig3); close(fig4);
display('Finished')
display('')
