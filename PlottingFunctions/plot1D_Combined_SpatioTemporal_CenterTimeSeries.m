function fig = plot1D_Combined_SpatioTemporal_CenterTimeSeries(BOLD_signal, ...
                deconvResponses, x, t, params, normalization, plot_what)
%% plot1D_Combined_SpatioTemporal_CenterTimeSeries.m
%
% Plots the 1D spatiotemporal responses and corresponding center responses 
% (at x = 0) in a signle figure.
%
% Inputs: BOLD_signal   : array of BOLD signal (x,t)
%                         size(BOLD_signal) = [length(distance), length(t)]
%         deconvResponses    : array of deconvolved responses (x,t for now)
%                         size(deconvResponses.{}) = [length(distance), length(t)]
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
% James Pang, University of Sydney, 2017

%%
% find index of t=0
t0 = dsearchn(t', 0);

% find index of x=0
x0 = dsearchn(x', 0);

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
    width = 0.25; height = 0.08; initial_x = 0.2; initial_y = 0.88;
    x_factor = 2; y_factor = 1.28;
    title_yloc = [2.5, 2, -3, 5, 5, 5, -2, -2, -2];

    fig = figure('Position', [200, 200, 400, 700], 'color','w');
    cmap = colormap_bluetored;
    
    for k=1:length(responses)-1
        if k==1
            data = real(BOLD_signal);
        else
            data = real(deconvResponses.(responses{k}));
        end
        subplot('Position', [initial_x initial_y-height*y_factor*(k-1) width*1.4 height]);
        imagesc(x*1e3, t(t0:end), data(:, t0:end).'); % note the transpose to make x the horizontal axis
        if k~=length(responses)-1
            set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], ...
                'ylim', [0, 20], 'xticklabel', {});
        else
            set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], ...
                'ylim', [0, 20]);
            xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
        end
        ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
        text(-12, title_yloc(k), titles{k}, 'fontsize', 15, 'rotation', 90)
    end
        
    for k=1:length(responses)-1
        if k==1
            data = real(BOLD_signal);
        else
            data = real(deconvResponses.(responses{k}));
        end
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));

        subplot('Position', [initial_x+width*x_factor initial_y-height*y_factor*(k-1) width height]);
        plot(t(t0:end), data(x0, t0:end), 'k-', 'linewidth', 2);
        if k~=length(responses)-1
            set(gca, 'fontSize', 13, 'xlim', [0, 20], 'ylim', [-clim_max, clim_max], 'xticklabel', {});
        else
            set(gca, 'fontSize', 13, 'xlim', [0, 20], 'ylim', [-clim_max, clim_max]);
            xlabel('$t$ (s)', 'fontsize',15, 'interpreter', 'latex')
        end
%             ylabel('response', 'fontsize',15)
    end
    
elseif strcmpi(plot_what, 'all_no_BOLD')
    width = 0.25; height = 0.087; initial_x = 0.2; initial_y = 0.86;
    x_factor = 2; y_factor = 1.28;
    title_yloc = [2.5, 2, -3, 5, 5, 5, -2, -2, -2];

    fig = figure('Position', [200, 200, 400, 600], 'color','w');
    cmap = colormap_bluetored;

    for k=2:length(responses)-1
        data = real(deconvResponses.(responses{k}));

        subplot('Position', [initial_x initial_y-height*y_factor*(k-2) width*1.4 height]);
        imagesc(x*1e3, t(t0:end), data(:, t0:end).'); % note the transpose to make x the horizontal axis
        if k~=length(responses)-1
            set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], ...
                'ylim', [0, 20], 'xticklabel', {});
        else
            set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], ...
                'ylim', [0, 20]);
            xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
        end
        ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
        text(-12, title_yloc(k), titles{k}, 'fontsize', 15, 'rotation', 90)
    end

    for k=2:length(responses)-1
        data = real(deconvResponses.(responses{k}));
        clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
        clim_max = max(abs(clim));

        subplot('Position', [initial_x+width*x_factor initial_y-height*y_factor*(k-2) width height]);
        plot(t(t0:end), data(x0, t0:end), 'k-', 'linewidth', 2);
        if k~=length(responses)-1
            set(gca, 'fontSize', 13, 'xlim', [0, 20], 'ylim', [-clim_max, clim_max], 'xticklabel', {});
        else
            set(gca, 'fontSize', 13, 'xlim', [0, 20], 'ylim', [-clim_max, clim_max]);
            xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
        end
%             ylabel('response', 'fontsize',15)
    end   
    
else
    width = 0.25; height = 0.5; initial_x = 0.2; initial_y = 0.25;
    x_factor = 2; 
    title_yloc = [4.5, 4, 1, 6, 6, 6, 1, 1, 1];
    
    fig = figure('Position', [200, 200, 400, 150], 'color','w');
    cmap = colormap_bluetored;
    
    response_index = find(strcmpi(responses, plot_what));
    if response_index==1
        data = real(BOLD_signal);
    else
        data = real(deconvResponses.(responses{response_index}));
    end

    subplot('Position', [initial_x initial_y width*1.4 height]);
    imagesc(x*1e3, t(t0:end), data(:, t0:end).'); % note the transpose to make x the horizontal axis
    set(gca, 'ydir', 'normal', 'fontSize', 13, 'xlim', [-5, 5], 'ylim', [0, 20]);
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    xlabel('$x$ (mm)','fontsize',15,'interpreter', 'latex')
    clim = [min(min(data(:,t0:end))), max(max(data(:,t0:end)))];
    clim_max = max(abs(clim));
    caxis([-clim_max, clim_max]);
    colormap(cmap)
    colorbar
    text(-12, title_yloc(response_index), titles{response_index}, ...
        'fontsize', 15, 'rotation', 90)


    subplot('Position', [initial_x+width*x_factor initial_y width height]);
    plot(t(t0:end), data(x0, t0:end), 'k-', 'linewidth', 2);
    set(gca, 'fontSize', 13, 'xlim', [0, 20], 'ylim', [-clim_max, clim_max]);
    xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
%         ylabel('response', 'fontsize',15)
end
