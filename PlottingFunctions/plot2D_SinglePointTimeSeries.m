function fig = plot2D_SinglePointTimeSeries(BOLD_signal, deconvResponses, x, y, t, ...
                    loc_interest, params, normalization, plot_what)
%% plot2D_SinglePointTimeSeries.m
%
% Plots the time series of responses at a particular location of interest
% in x,y space.

% Inputs: BOLD_signal   : array of BOLD signal (x,y,t)
%                         size(BOLD_signal) = [length(y), length(x), length(t)]
%         deconvResponses    : array of deconvolved responses (x,y,t)
%                         size(deconvResponses.{}) = [length(y), length(x), length(t)]
%                         Possible fields are reconvBOLD, neural, neuroglial, 
%                         CBF, CBV, dHb, Wmode, Lmode, and Dmode.
%         x             : vector of distance along x
%         y             : vector of distance along y
%         t             : vector of time 
%         loc_interest  : location of time series of interest in x,y space
%                         loc_interest(1) = point in x in mm
%                         loc_interest(2) = point in y in mm
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
% James Pang, University of Sydney, 2016

%%
% find index of location of interest
xp = dsearchn(x'*1e3, loc_interest(1));
yp = dsearchn(y'*1e3, loc_interest(2));

% find index of t=0
t0 = dsearchn(t', 0);

% find index of t=0.25
if BOLD_signal(yp, xp, t0)==0 && abs(diff(BOLD_signal(yp, xp, t0:t0+1))) > 0.05
    tstart = t0 + 1;
else
    tstart = t0;
end

% if true, normalize with respect to baseline value or maximum at t>0
if normalization
    BOLD_signal = BOLD_signal/max(max(max(real(BOLD_signal(:,:,t0:end)))));
    deconvResponses.reconvBOLD = deconvResponses.reconvBOLD/max(max(max(real(deconvResponses.reconvBOLD(:,:,t0:end)))));
    deconvResponses.neural = deconvResponses.neural/max(max(max(real(deconvResponses.neural(:,:,t0:end)))));
    deconvResponses.neuroglial = deconvResponses.neuroglial/max(max(max(real(deconvResponses.neuroglial(:,:,t0:end)))));
    deconvResponses.CBF = deconvResponses.CBF/params.F_0;
    deconvResponses.CBV = deconvResponses.CBV/params.Xi_0;
    deconvResponses.dHb = deconvResponses.dHb/params.Q_0;
    
    total_mode = real(deconvResponses.Wmode) + real(deconvResponses.Lmode) + ...
                 real(deconvResponses.Dmode);
    deconvResponses.Wmode = deconvResponses.Wmode/max(max(max(total_mode(:,:,t0:end))));
    deconvResponses.Lmode = deconvResponses.Lmode/max(max(max(total_mode(:,:,t0:end))));
    deconvResponses.Dmode = deconvResponses.Dmode/max(max(max(total_mode(:,:,t0:end))));
end

titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode', 'reconv BOLD'};
responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode', 'reconvBOLD'};

if strcmpi(plot_what, 'all_w_BOLD')
    width = 0.188; height = 0.2; initial_x = 0.055; initial_y = 0.72;
    x_factor = 1.28; y_factor = 1.55;

    fig = figure('Position', [200, 200, 800, 500]);
    
    data = squeeze(real(BOLD_signal(yp,xp,tstart:end)));
    subplot(3,4,2, 'Parent', fig, 'Position', [0.43 initial_y width height]);
    plot(t(tstart:end), data, 'k-', 'Linewidth', 2);  
    title(titles{1}, 'fontsize', 15);
    set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
        'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%     xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')

    for sub = 5:8
        data = squeeze(real(deconvResponses.(responses{sub-3})(yp,xp,tstart:end)));
        subplot(3,4,sub, 'Parent', fig, 'Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height]);
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub-3}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%         xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end

    for sub = 9:12
        data = squeeze(real(deconvResponses.(responses{sub-3})(yp,xp,tstart:end)));
        subplot(3,4,sub, 'Parent', fig, 'Position', [initial_x+width*x_factor*(sub-9) initial_y-height*y_factor*2 width height])
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
        data = squeeze(real(deconvResponses.(responses{sub+1})(yp,xp,tstart:end)));
        subplot(2,4,sub, 'Parent', fig, 'Position', [initial_x+width*x_factor*(sub-1) initial_y width height]);
        plot(t(tstart:end), data, 'k-', 'Linewidth', 2);
        title(titles{sub+1}, 'fontsize',15);
        set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
            'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
%         xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    end
    
    for sub = 5:8
        data = squeeze(real(deconvResponses.(responses{sub+1})(yp,xp,tstart:end)));
        subplot(2,4,sub, 'Parent', fig, 'Position', [initial_x+width*x_factor*(sub-5) initial_y-height*y_factor width height])
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
        data = squeeze(real(BOLD_signal(yp,xp,tstart:end)));
    else
        data = squeeze(real(deconvResponses.(responses{response_index})(yp,xp,tstart:end)));
    end
    
    plot(t(tstart:end), data, 'k-', 'Linewidth', 2);  
    set(gca, 'fontSize', 13, 'xlim', [0, 20], 'xtick', 0:5:20, ...
        'ylim', [min(data(:))*1.1, max(data(:))*1.1]);
    xlabel('$t$ (s)','fontsize',15,'interpreter', 'latex')
    ylabel(titles{response_index},'fontsize',15) 
end
