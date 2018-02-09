function fig = plot1D_SpatioTemporal(BOLD_signal, deconvResponses, x, t, params, ...
                             normalization, plot_what)
%% plot1D_SpatioTemporal.m
%
% Plots the 1D spatiotemporal responses as contour plots
%
% Inputs: BOLD_signal   : array of BOLD signal (x,t)
%                         size(BOLD_signal) = [length(distance), length(t)]
%         deconvResponses    : array of deconvolved responses (x,t)
%                         size(responses.{}) = [length(distance), length(t)]
%                         Possible fields are reconvBOLD, neural, neuroglial, 
%                         CBF, CBV, dHb, Wmode, Lmode, and Dmode.
%         x             : vector of distance
%         t             : vector of time 
%         params        : instance of the class loadParameters of the 
%                         toolbox
%         normalization : 1 or 0
%                         Choose 1 if you want the responses normalized, 
%                         0 otherwise.
%         plot_what     : string saying which response(s) to plot.
%                         Possible inputs are 'BOLD', 'reconvBOLD', 'neural', 
%                         'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 
%                         'Dmode', 'all_w_BOLD', 'all_no_BOLD'.
%
% Output: fig           : figure handle of the resulting plot
%
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%
% find index of t=0
t0 = dsearchn(t', 0);

% if true, normalize with respect to baseline value or maximum at t>0
if normalization
    BOLD_signal = BOLD_signal/max(max(real(BOLD_signal(:,t0:end))));
    deconvResponses.reconvBOLD = deconvResponses.reconvBOLD/max(max(max(real(deconvResponses.reconvBOLD(:,:,t0:end)))));
    deconvResponses.neural = deconvResponses.neural/max(max(max(real(deconvResponses.neural(:,:,t0:end)))));
    deconvResponses.neuroglial = deconvResponses.neuroglial/max(max(real(deconvResponses.neuroglial(:,t0:end))));
    deconvResponses.CBF = deconvResponses.CBF/params.F_0;
    deconvResponses.CBV = deconvResponses.CBV/params.Xi_0;
    deconvResponses.dHb = deconvResponses.dHb/params.Q_0;
    
    total_mode = real(deconvResponses.Wmode) + real(deconvResponses.Lmode) + ...
                 real(deconvResponses.Dmode);
    deconvResponses.Wmode = deconvResponses.Wmode/max(max(total_mode(:,t0:end)));
    deconvResponses.Lmode = deconvResponses.Lmode/max(max(total_mode(:,t0:end)));
    deconvResponses.Dmode = deconvResponses.Dmode/max(max(total_mode(:,t0:end)));
end

titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode', 'reconv BOLD'};
responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode', 'reconvBOLD'};
         
if strcmpi(plot_what, 'all_w_BOLD')
    width = 0.2; height = 0.2; initial_x = 0.055; initial_y = 0.72;
    x_factor = 1.2; y_factor = 1.55;

    fig = figure('Position', [200, 200, 800, 500]);
    cmap = colormap_bluetored;
    
    data = real(BOLD_signal);
    subplot('Position', [0.43 initial_y width height]);
    imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
    title(titles{1}, 'fontsize', 15);
    set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
    clim_max = max(abs(clim));
    caxis([-clim_max, clim_max]);
    colormap(cmap)
    colorbar
    
    for sub = 5:8
        data = real(deconvResponses.(responses{sub-3}));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height]);
        imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
        title(titles{sub-3}, 'fontsize',15);
        set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
        if sub==5
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end

    for sub = 9:12
        data = real(deconvResponses.(responses{sub-3}));
        subplot('Position', [initial_x+width*x_factor*(sub-9) initial_y-height*y_factor*2 width height])
        imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
        title(titles{sub-3},'fontsize',15);
        set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
        xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
        if sub==9
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end
    
elseif strcmpi(plot_what, 'all_no_BOLD')
    width = 0.2; height = 0.3; initial_x = 0.055; initial_y = 0.58;
    x_factor = 1.2; y_factor = 1.5;

    fig = figure('Position', [200, 200, 800, 400]);
    cmap = colormap_bluetored;
    
    for sub = 1:4
        data = real(deconvResponses.(responses{sub+1}));
        subplot('Position', [initial_x+width*x_factor*(sub-1) initial_y width height]);
        imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
        title(titles{sub+1}, 'fontsize',15);
        set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
        if sub==1
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end

    for sub = 5:8
        data = real(deconvResponses.(responses{sub+1}));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height])
        imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
        title(titles{sub+1},'fontsize',15);
        set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
        xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
        if sub==5
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end
    
else
    fig = figure('Position', [200, 200, 400, 300]);
    cmap = colormap_bluetored;
    
    response_index = find(strcmpi(responses, plot_what));
    if response_index==1
        data = real(BOLD_signal);
    else
        data = real(deconvResponses.(responses{response_index}));
    end
    
    imagesc(x*1e3, t, data.'); % note the transpose to make x the horizontal axis
    title(titles{response_index}, 'fontsize',15);
    set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
    xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
    clim_max = max(abs(clim));
    caxis([-clim_max, clim_max]);
    colormap(cmap)
    colorbar
end        
