function fig = plot_V1TimeEvolution(responses_1D, eccentricity, t, params, ...
                                    normalization, clim_factor, plot_what)
%% plot_V1TimeEvolution.m
%
% Plots the average evolution of the responses in V1 as a function
% of time and eccentricity
%
% Inputs: responses_1D  : array of average 1D responses in V1 BOLD (eccentricity, t)
%                         size(responses_1D.{}) = [length(eccentricity), length(t)]
%                         Possible fields are BOLD, reconvBOLD, neural, 
%                         neuroglial, CBF, CBV, dHb, Wmode, Lmode, and Dmode.
%         eccentricity  : vector of eccentricity values
%         t             : vector of time 
%         params        : instance of the class loadParameters of the 
%                         toolbox
%         normalization : 1 or 0
%                         choose 1 if you want the responses normalized
%         clim_factor   : dividing the maximum value of clim by this value
%         plot_what     : string saying which response to plot.
%                         Possible inputs are 'BOLD', 'reconvBOLD', 'neural', 
%                         'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 
%                         'Dmode', 'all_w_BOLD', 'all_no_BOLD'
%
% Output: fig           : figure handle of the resulting plot
%
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% if true, normalize with respect to baseline value or maximum 
if normalization
    responses_1D.BOLD = responses_1D.BOLD/max(max(real(responses_1D.BOLD)));
    responses_1D.reconvBOLD = responses_1D.reconvBOLD/max(max(real(responses_1D.reconvBOLD)));
    responses_1D.neural = responses_1D.neural/max(max(real(responses_1D.neural)));
    responses_1D.neuroglial = responses_1D.neuroglial/max(max(real(responses_1D.neuroglial)));
    responses_1D.CBF = responses_1D.CBF/params.F_0;
    responses_1D.CBV = responses_1D.CBV/params.Xi_0;
    responses_1D.dHb = responses_1D.dHb/params.Q_0;
    
    total_mode = real(responses_1D.Wmode) + real(responses_1D.Lmode) + ...
                 real(responses_1D.Dmode);
    responses_1D.Wmode = responses_1D.Wmode/max(total_mode(:));
    responses_1D.Lmode = responses_1D.Lmode/max(total_mode(:));
    responses_1D.Dmode = responses_1D.Dmode/max(total_mode(:));
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
    
    data = real(responses_1D.(responses{1}));
    subplot('Position', [0.43 initial_y width height]);
    contourf(eccentricity, t, data, 10, 'EdgeColor','none'); 
    title(titles{1}, 'fontsize', 15);
    set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
        'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
        'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    clim = [min(data(:)), max(data(:))];
    clim_max = max(abs(clim))/clim_factor;
    caxis([-clim_max, clim_max]);
    colormap(cmap)
    colorbar

    for sub = 5:8
        data = real(responses_1D.(responses{sub-3}));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height]);
        contourf(eccentricity, t, data, 10, 'EdgeColor','none');
        title(titles{sub-3}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
            'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
            'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
        if sub==5
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(data(:)), max(data(:))];
        clim_max = max(abs(clim))/clim_factor;
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end

    for sub = 9:12
        data = real(responses_1D.(responses{sub-3}));
        subplot('Position', [initial_x+width*x_factor*(sub-9) initial_y-height*y_factor*2 width height])
        contourf(eccentricity, t, data, 10, 'EdgeColor','none');
        title(titles{sub-3},'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
            'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
            'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
        xlabel('eccentricity ($^\circ$)','fontsize',15,'interpreter', 'latex')
        if sub==9
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(data(:)), max(data(:))];
        clim_max = max(abs(clim))/clim_factor;
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
        data = real(responses_1D.(responses{sub+1}));
        subplot('Position', [initial_x+width*x_factor*(sub-1) initial_y width height]);
        contourf(eccentricity, t, data, 10, 'EdgeColor','none');
        title(titles{sub+1}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
            'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
            'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
        if sub==1
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(data(:)), max(data(:))];
        clim_max = max(abs(clim))/clim_factor;
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end

    for sub = 5:8
        data = real(responses_1D.(responses{sub+1}));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height])
        contourf(eccentricity, t, data, 10, 'EdgeColor','none');
        title(titles{sub+1},'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
            'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
            'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
        xlabel('eccentricity ($^\circ$)','fontsize',15,'interpreter', 'latex')
        if sub==5
            ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
        end
        clim = [min(data(:)), max(data(:))];
        clim_max = max(abs(clim))/clim_factor;
        caxis([-clim_max, clim_max]);
        colormap(cmap)
        colorbar
    end
    
else
    fig = figure('Position', [200, 200, 400, 300]);
    cmap = colormap_bluetored;
    
    response_index = find(strcmpi(responses, plot_what));
    data = real(responses_1D.(responses{response_index}));
    
    contourf(eccentricity, t, data, 10, 'EdgeColor','none');
    title(titles{response_index}, 'fontsize',15);
    set(gca, 'fontSize', 13, 'xlim', [eccentricity(1), eccentricity(end)], ...
            'ylim', [t(1), t(end)], 'xtick', [eccentricity(1), 1:4, eccentricity(end)], ...
            'xticklabel', [0, 1, 2, 3, 4, eccentricity(end)]);
    xlabel('eccentricity ($^\circ$)','fontsize',15,'interpreter', 'latex')
    ylabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    clim = [min(data(:)), max(data(:))];
    clim_max = max(abs(clim))/clim_factor;
    caxis([-clim_max, clim_max]);
    colormap(cmap)
    colorbar
end
