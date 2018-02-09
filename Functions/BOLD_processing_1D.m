function [BOLD_processed, x, t] = BOLD_processing_1D(BOLD_signal, ...
                                    x_experiment, t_experiment, params) 
%% BOLD_processing_1D.m
%
% Performs processing works on the 1D BOLD signal. This function
% interpolates the BOLD signal to a desired resolution based on params.Nkx
% and params.Nw. It then pads zeros to the temporal dimension of the
% interpolated BOLD signal to make it centered at t=0 and make BOLD(at t<0)
% = 0;
%
% Inputs: BOLD_signal       : array of 1D BOLD signal (x,t)
%                             size(BOLD_signal) = [length(x_experiment), 
%                                                  length(t_experiment)] 
%         x_experiment      : vector of distance to get BOLD_signal 
%         t_experiment      : vector of time to get BOLD_signal 
%         params            : instance of the class loadParameters of the 
%                             toolbox
%
% Output: BOLD_processed    : array of interpolated and zero-padded BOLD signal.
%                             size(BOLD_signal) = [params.Nkx, params.Nw]
%         x                 : vector of new distance x =[-x_end,...,0,...x_end]
%         t                 : vector of new time t =[-t_end,...,0,...t_end]
% 
% Example:
% >> params = loadParameters;
% >> load BOLD_signal.mat   % assuming the data is stored in this mat file
% >> x_experiment = linspace(-5,5,256)*1e-3;  % in mm
% >> t_experiment = linspace(0.1,20,256);     % in s
% >> [BOLD_processed, x, t] = BOLD_processing_1D(BOLD_signal, x_experiment, 
%                     t_experiment, params) % gives out the processed BOLD
%                                             and new time/distance vectors 
% 
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%
% interpolate experimental response to increase resolution

% if the start of t_experiment has an offset from zero, t_interp becomes of
% size (Nw-1)/2 and then add a zero in the middle of t
% if the start of t_experiment is zero, t_interp becomes of size (Nw+1)/2 
% and there is no need to add a zero in the middle of t 
if t_experiment(1)~=0
    t_interp = linspace(t_experiment(1), t_experiment(end), params.Nw/2);
    t = [-t_interp(end:-1:1), 0, t_interp];
else
    t_interp = linspace(t_experiment(1), t_experiment(end), params.Nw/2 + 1);
    t = [-t_interp(end:-1:2), t_interp];
end

% assuming that x_experiment is symmetric about zero and  
% x_experiment = -x_experiment(end)
x_interp = x_experiment(end)*(2*(0:params.Nkx)/params.Nkx - 1);

% creating matrices of time and distance
[t_interp_mat, x_interp_mat] = meshgrid(t_interp, x_interp);
[t_experiment_mat, x_experiment_mat] = meshgrid(t_experiment, x_experiment);

% 2D interpolation
BOLD_signal_interp = interp2(t_experiment_mat, x_experiment_mat, BOLD_signal, ...
                             t_interp_mat, x_interp_mat);

% since we increased the resolution of time from -t to t, we need to pad 
% the interpolated BOLD signal with zeros such that BOLD(t<=0) = 0                         
BOLD_signal_padded = padarray(BOLD_signal_interp, ...
                             [0, length(t)-length(t_interp)], 'pre');

% reassigning the resulting interpolated and zero-padded BOLD signal to the 
% final variable BOLD_processed                             
BOLD_processed = BOLD_signal_padded;

% Fourier transform requires that the x=0 and t=0 be in an offset to the
% right such that an N-sized vector F, where N is even, needs to have the 
% property F(1) = F(end-1) and F(N/2+1) corresponds to the center (either
% the x=0 or t=0)
x = x_interp(1:end-1);
t = t(1:end-1);
BOLD_processed = BOLD_processed(1:end-1, 1:end-1);
