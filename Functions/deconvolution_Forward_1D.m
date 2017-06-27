function deconvResponses = deconvolution_Forward_1D(neural_signal, ...
                                            distance, time, params) 
%% deconvolution_Forward_1D.m
%
% Deconvolves the 1D responses using the Forward method described in the paper.
% The forward method assumes that an estimated neural activity is
% available.
%
% Inputs: neural_signal : array of 1D neural activity (x,t)
%                         size(neural_signal) = [length(distance), length(time)] 
%         distance      : vector of distance to get neural_signal. 
%                         This needs to be symmetric with respect to x = 0 
%                         such that x(1) = -x(end)
%         time          : vector of time to get neural_signal. 
%                         This needs to be symmetric with respect to t = 0 
%                         such that t(1) = -t(end).
%         params        : instance of the class loadParameters of the 
%                         toolbox
%
% Output: deconvResponses    : structure containing the 1D deconvolved responses.
%                         Possible fields are reconvBOLD, neural, neuroglial, 
%                         CBF, CBV, dHb, Wmode, Lmode, and Dmode.
% 
% Example:
% >> params = loadParameters;
% >> load neural_signal.mat               % assuming neural data is stored in this mat file
% >> distance = linspace(-5,5,256)*1e-3;  % in mm
% >> time = linspace(-20,20,256);         % in s
% >> deconvResponses = deconvolution_Forward_1D(neural_signal, distance, time, params)
% >> deconvResponses.CBF                  % gives out the deconvolved 1D CBF
% 
% James Pang, University of Sydney, 2016

%%

% calculating computational parameters 
Nkx = length(distance);
Nw = length(time);
kxsamp = (1/mean(diff(distance)))*2*pi;
wsamp = (1/mean(diff(time)))*2*pi;

% creating vectors of frequencies kx and w
[kx, w] = generate_kw_1D(kxsamp, wsamp, Nkx, Nw);

% calculate 1D transfer function from phi
% T is a structure with possible fields: T_Yphi, T_zphi, T_Fphi, T_Xiphi, 
% T_Qphi, T_1phi, T_2phi, T_3phi, T_4phi, T_5phi, T_Wphi, T_Lphi, and T_Dphi.
T = calcTransFuncs_fromPhi_1D(kx, w, params);

% assigning neural signal in the structure and calculating its Fourier
% transform
deconvResponses.('neural') = neural_signal;
neural_signal_freq = coord2freq_1D(neural_signal, kx, w);

% deconvolve 1D responses using a forward calculation  
names = {'reconvBOLD', 'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 'Dmode', 'BOLD'};
variables = {'Y', 'z', 'F', 'Xi', 'Q', 'W', 'L', 'D', 'Y'};
for i=1:length(variables)
    trans_func = ['T.T_', variables{i}, 'phi'];
    deconvResponses.(names{i}) = real(freq2coord_1D(eval(trans_func).*neural_signal_freq, ...
                            kx, w));
end
