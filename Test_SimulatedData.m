%% Test_SimulatedData.m
%
% This test script shows how the deconvolution method works for a 1D 
% simulated BOLD data driven by a noisy spatiotemporal stationary Gaussian 
% neural activity
%
% James Pang, University of Sydney, 2016

%% Adding the paths of the sub-directories for direct access of files
% This is not necessary if entire Palimpsest toolbox is added to the Matlab 
% path via addpath(genpath('PalimpsestToolboxLocation')) where
% PalimpsestToolboxLocation is the location of the toolbox

% addpath('Data', 'Functions', 'PlottingFunctions')

%% Loading the default values of the model parameters 

params = loadParameters;

%% Initializing computational parameters

d_max = 15e-3;
t_max = 20;

params.Nkx = 2^9;       % increase resolution for 1D version
params.Nw = 2^11;       % increase resolution for 1D version

% distance and time, in 1D
% note that Fourier transform requires that the vectors be defined this way 
% and not via using linspace 
x = d_max*(2*(0:params.Nkx-1)/params.Nkx - 1);
t = t_max*(2*(0:params.Nw-1)/params.Nw - 1);

kxsamp = (1/mean(diff(x)))*2*pi;
wsamp = (1/mean(diff(t)))*2*pi;

[kx, w] = generate_kw_1D(kxsamp, wsamp, params.Nkx, params.Nw);

%% Generating simulated neural activity and BOLD 
[tt, xx] = meshgrid(t, x);
sigma_t = 1;
sigma_x = 1e-3;
t_0 = 2;                % fixed activity time offset

UseSampleData = 1;

if UseSampleData
    load('Data/Simulated/simulated_1Ddata.mat', 'simulated_phi_noisy', 'simulated_BOLD_noisy');
else
    noise_term = 0.1*randn(size(tt));
    simulated_phi = exp(-(tt-t_0).^2/sigma_t^2).*exp(-xx.^2/sigma_x^2);
    simulated_phi_noisy = simulated_phi + noise_term;
    simulated_phi_noisy(tt<=0) = 0;

    T = calcTransFuncs_fromPhi_1D(kx, w, params);

    ds = mean(diff(t))*mean(diff(x))*1e3;
    stHRF = real(freq2coord_1D(T.T_Yphi, kx, w));
    simulated_BOLD = ds*real(convnfft(stHRF, simulated_phi, 'same'));
    simulated_BOLD_noisy = ds*real(convnfft(stHRF, simulated_phi_noisy, 'same'));

    save('Data/Simulated/simulated_1Ddata.mat', 'simulated_phi_noisy', 'simulated_BOLD_noisy')
end

%% Deconvolution of responses 

% changing the NSR term
% params.noise = [300, 0.1];         % cut-off frequencies
params.noise = 0.5;                % constant

deconvResponses_Forward_1D = deconvolution_Forward_1D(simulated_phi_noisy, ...
                                            x, t, params); 
deconvResponses_HybridWiener_1D = deconvolution_HybridWiener_1D(simulated_BOLD_noisy, ...
                                            x, t, params);

%% Plotting the spatiotemporal results in a contour 

titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode'};
responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode'};

width = 0.3; height = 0.072; initial_x = 0.26; initial_y = 0.88;
x_factor = 1.2; y_factor = 1.35;

% fix limits
cmin = [0, 0, 0, -0.4, -50, -0.025, -0.16, -0.15, -0.18];
cmax = [0.64, 0.6, 0.6, 1.6, 232, 0.005, 0.15, 0.7, 0.01];

% find index of t=0
t0 = dsearchn(t', 0);

fig1 = figure('Position', [200, 200, 400, 1000], 'color','w');
cmap_new = colormap_bluetored;

data = simulated_phi_noisy;
subplot('Position', [0.2 initial_y width height])
imagesc(x*1e3, t, data.');
title('simulated neural', 'fontsize', 15);
set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [min(x)*1e3, max(x)*1e3], ...
    'ylim', [0, 15], 'ticklength', 3*get(gca,'ticklength'));
% xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
clim = [cmin(2), cmax(2)];
clim_max = max(abs(clim));
caxis([-clim_max, clim_max]);
colormap(cmap_new)
colorbar
annotation('textbox', get(gca, 'Position')+[-0.13,0.045,0,0], 'String', '(a)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')

data = simulated_BOLD_noisy;
subplot('Position', [0.2+width*x_factor*1.2 initial_y width height]);
imagesc(x*1e3, t, data.');
title('simulated BOLD', 'fontsize', 15);
set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [min(x)*1e3, max(x)*1e3], ...
    'ylim', [0, 15], 'ticklength', 3*get(gca,'ticklength'));
% xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
clim = [cmin(1), cmax(1)];
clim_max = max(abs(clim));
caxis([-clim_max, clim_max]);
colormap(cmap_new)
colorbar
annotation('textbox', get(gca, 'Position')+[-0.13,0.045,0,0], 'String', '(b)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')

for sub=8:15
    data = real(deconvResponses_Forward_1D.(responses{sub-6}));
    subplot('Position', [initial_x initial_y-height*y_factor*(sub-6.5) width height])
    imagesc(x*1e3, t, data.');

    set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [min(x)*1e3, max(x)*1e3], ...
        'ylim', [0, 15], 'ticklength', 3*get(gca,'ticklength'));
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    if sub==8
        annotation('textbox', get(gca, 'Position')+[0.04,0.04,0,0], 'String', '(c)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    end
    if sub==8 || sub==10 || sub==11 || sub==12
        text(min(x)*1e3*2.7, 3, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    else
        text(min(x)*1e3*2.7, -4, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    end
    if sub~=15
        set(gca, 'xticklabel', {})
    else
        xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    end
    clim = [cmin(sub-6), cmax(sub-6)];
    clim_max = max(abs(clim));
    caxis([-clim_max, clim_max]);
    colormap(cmap_new)
    colorbar
end

for sub=16:23 
    data = real(deconvResponses_HybridWiener_1D.(responses{sub-14}));

    subplot('Position', [initial_x+width*x_factor initial_y-height*y_factor*(sub-14.5) width height])
    imagesc(x*1e3, t, data.');

    set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [min(x)*1e3, max(x)*1e3], ...
        'ylim', [0, 15], 'ticklength', 3*get(gca,'ticklength'));
    if sub==16
        annotation('textbox', get(gca, 'Position')+[0.04,0.04,0,0], 'String', '(d)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    end
    if sub~=23
        set(gca, 'xticklabel', {}, 'yticklabel', {})
    else
        set(gca, 'yticklabel', {})
        xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    end
    clim = [cmin(sub-14), cmax(sub-14)];
    clim_max = max(abs(clim));
    caxis([-clim_max, clim_max]);
    colormap(cmap_new)
    colorbar
end

set(fig1, 'PaperPositionMode','auto')     %# WYSIWYG
print(fig1, '-painters', '-depsc', 'Figures/SimulatedResults_SpatioTemporal.eps')

%% Plotting the time series of center response (at x = 0) results 

titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode'};
responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode'};

width = 0.25; height = 0.069; initial_x = 0.24; initial_y = 0.88;
x_factor = 1.5; y_factor = 1.4;

% find index of x=0
center = dsearchn(x', 0);

% fix limits
ymin = [-0.2, -0.05, -0.05, -0.28, -70.3, -0.025, -0.17, -0.16, -0.19];
ymax = [0.85, 1.2, 1.2, 1.6, 235, 0.007, 0.032, 0.77, 0.013];

fig2 = figure('Position', [200, 200, 400, 1000]);

data = real(simulated_phi_noisy(center,:));
subplot('Position', [0.21 initial_y width height])
plot(t, data, 'k-', 'Linewidth', 2);
title('simulated neural', 'fontsize', 15);
set(gca, 'fontSize', 13, 'xlim', [0, 15], 'xtick', 0:5:20, ...
        'ylim', [ymin(2), ymax(2)], 'ticklength', 3*get(gca,'ticklength'));
% xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
annotation('textbox', get(gca, 'Position')+[-0.1,0.045,0,0], 'String', '(a)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')

data = real(simulated_BOLD_noisy(center,:));
subplot('Position', [0.21+width*x_factor*1.1 initial_y width height]);
plot(t, data, 'k-', 'Linewidth', 2);
title('simulated BOLD', 'fontsize', 15);
set(gca, 'fontSize', 13, 'xlim', [0, 15], 'xtick', 0:5:20, ...
        'ylim', [ymin(1), ymax(1)], 'ytick', [0,0.4,0.8],'ticklength', 3*get(gca,'ticklength'));
% xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
annotation('textbox', get(gca, 'Position')+[-0.1,0.045,0,0], 'String', '(b)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')

for sub=8:15
    data = real(deconvResponses_Forward_1D.(responses{sub-6})(center,:));
    subplot('Position', [initial_x initial_y-height*y_factor*(sub-6.5) width height])
    plot(t, data, 'k-', 'Linewidth', 2)
    set(gca, 'fontSize', 13, 'xlim', [0, 15], 'xtick', 0:5:20, ...
        'ylim', [ymin(sub-6), ymax(sub-6)], 'ticklength', 3*get(gca,'ticklength'));
    if sub==8
        annotation('textbox', get(gca, 'Position')+[0.065,0.04,0,0], 'String', '(c)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    end
    if sub==9
        text(-10, ymin(sub-6)*10, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    elseif sub==8 || sub==10 || sub==11 || sub==12
        text(-10, ymin(sub-6)*0.8, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    elseif sub==13 
        text(-10, ymin(sub-6)*1.35, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    elseif sub==14 
        text(-10, ymin(sub-6)*2, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    elseif sub==15
        text(-10, ymin(sub-6)*1.3, titles{sub-6}, 'fontsize', 15, 'rotation', 90)
    end
    if sub~=15
        set(gca, 'xticklabel', {})
    else 
        xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end
end
% text(-13, 0.4, 'Deconvolved Responses', 'fontsize', 15, 'rotation', 90, 'fontweight', 'b')

for sub=16:23
    data = real(deconvResponses_HybridWiener_1D.(responses{sub-14})(center,:));
    subplot('Position', [initial_x+width*x_factor initial_y-height*y_factor*(sub-14.5) width height])
    plot(t, data, 'k-', 'Linewidth', 2)
    set(gca, 'fontSize', 13, 'xlim', [0, 15], 'xtick', 0:5:20, ...
        'ylim', [ymin(sub-14), ymax(sub-14)], 'ticklength', 3*get(gca,'ticklength'));
    if sub==16
        annotation('textbox', get(gca, 'Position')+[0.065,0.04,0,0], 'String', '(d)', 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    end
    if sub~=23
        set(gca, 'xticklabel', {})
    else
        xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end
end

set(fig2, 'PaperPositionMode','auto')     %# WYSIWYG
print(fig2, '-painters', '-depsc', 'Figures/SimulatedResults_CenterTimeSeries.eps')
