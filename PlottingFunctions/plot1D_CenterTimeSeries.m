function fig = plot1D_CenterTimeSeries(BOLD_signal, deconvResponses, x, t, params, ...
                             normalization, plot_what)
%% plot1D_CenterTimeSeries.m
%
% Plots the time series of the 1D center responses (at x = 0).
%
% Inputs: BOLD_signal   : array of 1D BOLD signal (x,t)
%                         size(BOLD_signal) = [length(distance), length(t)]
%         deconvResponses    : array of 1D deconvolved responses (x,t)
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
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%
% find index of x=0
x0 = dsearchn(x', 0);

% find index of t=0
t0 = dsearchn(t', 0);

% find index of t=0.25
if BOLD_signal(x0, t0)==0 && abs(diff(BOLD_signal(x0, t0:t0+1))) > 0.05
    tstart = t0 + 1;
else
    tstart = t0;
end

% if true, normalize with respect to baseline value or maximum at t>0
if normalization
    BOLD_signal = BOLD_signal/max(max(real(BOLD_signal(:,tstart:end))));
    deconvResponses.reconvBOLD = deconvResponses.reconvBOLD/max(max(max(real(deconvResponses.reconvBOLD(:,:,t0:end)))));
    deconvResponses.neural = deconvResponses.neural/max(max(max(real(deconvResponses.neural(:,:,t0:end)))));
    deconvResponses.neuroglial = deconvResponses.neuroglial/max(max(real(deconvResponses.neuroglial(:,tstart:end))));
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
    width = 0.188; height = 0.2; initial_x = 0.055; initial_y = 0.72;
    x_factor = 1.28; y_factor = 1.55;

    fig = figure('Position', [200, 200, 800, 500]);
    
    data = real(BOLD_signal(x0,tstart:end));
    subplot('Position', [0.43 initial_y width height]);
    plot(t(tstart:end), data, 'k-', 'Linewidth', 2);  
    title(titles{1}, 'fontsize', 15);
    set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
        'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%     xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')

    for sub = 5:8
        data = real(deconvResponses.(responses{sub-3})(x0,tstart:end));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height]);
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub-3}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%         xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end

    for sub = 9:12
        data = real(deconvResponses.(responses{sub-3})(x0,tstart:end));
        subplot('Position', [initial_x+width*x_factor*(sub-9) initial_y-height*y_factor*2 width height])
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub-3},'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
        xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end
    
elseif strcmpi(plot_what, 'all_no_BOLD')
    width = 0.188; height = 0.28; initial_x = 0.06; initial_y = 0.58;
    x_factor = 1.28; y_factor = 1.55;

    fig = figure('Position', [200, 200, 800, 400]);
    
    for sub = 1:4
        data = real(deconvResponses.(responses{sub+1})(x0,tstart:end));
        subplot('Position', [initial_x+width*x_factor*(sub-1) initial_y width height]);
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub+1}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%         xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end

    for sub = 5:8
        data = real(deconvResponses.(responses{sub+1})(x0,tstart:end));
        subplot('Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height])
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub+1},'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
        xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end
        
else
    fig = figure('Position', [200, 200, 400, 300]);
    
    response_index = find(strcmpi(responses, plot_what));
    if response_index==1
        data = real(BOLD_signal(x0,tstart:end));
    else
        data = real(deconvResponses.(responses{response_index})(x0,tstart:end));
    end
    
    plot(t(tstart:end), data, 'k-', 'Linewidth', 2);  
    set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
        'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
    xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    ylabel(titles{response_index},'fontsize',15) 
end
