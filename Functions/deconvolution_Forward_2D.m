function deconvResponses = deconvolution_Forward_2D(neural_signal, ...
                                            distancex, distancey, time, params) 
%% deconvolution_Forward_2D.m
%
% Deconvolves the 2D responses using the Forward method described in the paper.
% The forward method assumes that an estimated neural activity is
% available.
%
% Inputs: BOLD_signal   : array of 2D BOLD signal (x,y,t)
%                         size(BOLD_signal) = [length(distancey), 
%                                          length(distancex), length(time)] 
%         neural_signal : array of 2D neural activity (x,y,t)
%                         size(neural_signal) = [length(distancey), 
%                                          length(distancex), length(time)] 
%         distancex     : vector of distance along x to get the neural_signal. 
%                         This needs to be symmetric with respect to x = 0 
%                         such that x(1) = -x(end)
%         distancey     : vector of distance along y to get the neural_signal. 
%                         This needs to be symmetric with respect to y = 0 
%                         such that y(1) = -y(end)
%         time          : vector of time to get neural_signal. 
%                         This needs to be symmetric with respect to t = 0 
%                         such that t(1) = -t(end).
%         params        : instance of the class loadParameters of the 
%                         toolbox
%
% Output: deconvResponses    : structure containing the 2D deconvolved responses.
%                         Possible fields are reconvBOLD, neural, neuroglial, 
%                         CBF, CBV, dHb, Wmode, Lmode, and Dmode.
% 
% Example:
% >> params = loadParameters;
% >> load neural_signal.mat   % assuming neural data is stored in this mat file
% >> distancex = linspace(-5,5,256)*1e-3;  % in mm
% >> distancey = linspace(-5,5,256)*1e-3;  % in mm
% >> time = linspace(-20,20,256);          % in s
% >> deconvResponses = deconvolution_Forward_2D(neural_signal,  
%                                       distancex, distancey, time, params)
% >> deconvResponses.CBF                   % gives out the deconvolved 2D CBF
% 
% James Pang, University of Sydney, 2016

%%

% calculating computational parameters 
Nkx = length(distancex);
Nky = length(distancey);
Nw = length(time);

kxsamp = (1/mean(diff(distancex)))*2*pi;
kysamp = (1/mean(diff(distancey)))*2*pi;
wsamp = (1/mean(diff(time)))*2*pi;

% creating vectors of frequencies kx, ky, and w
[kx, ky, w] = generate_kw_2D(kxsamp, kysamp, wsamp, Nkx, Nky, Nw);

% calculate 2D transfer function from Y
% T is a structure with possible fields: T_Yphi, T_zphi, T_Fphi, T_Xiphi, 
% T_Qphi, T_1phi, T_2phi, T_3phi, T_4phi, T_5phi, T_Wphi, T_Lphi, and T_Dphi.
T = calcTransFuncs_fromPhi_2D(kx, ky, w, params);

% assigning neural signal in the structure and calculating its Fourier
% transform
deconvResponses.('neural') = neural_signal;
neural_signal_freq = coord2freq_2D(neural_signal, kx, ky, w);

% deconvolve 2D responses using a forward calculation 
names = {'reconvBOLD', 'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 'Dmode'};
variables = {'Y', 'z', 'F', 'Xi', 'Q', 'W', 'L', 'D'};
for i=1:length(variables)
    trans_func = ['T.T_', variables{i}, 'phi'];
    deconvResponses.(names{i}) = real(freq2coord_2D(eval(trans_func).*neural_signal_freq, ...
                            kx, ky, w));
end
