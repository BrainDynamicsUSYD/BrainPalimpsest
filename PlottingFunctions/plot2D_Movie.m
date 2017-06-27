function fig = plot2D_Movie(BOLD_signal, deconvResponses, x, y, t, ...
                    zoom_lim, params, normalization, clim_factor, plot_what,...
                    movie_filename, frame_rate, v1_boundary)
%% plot2D_Movie.m
%
% Plots the evolution of the 2D spatiotemporal responses as a movie
%
% Inputs: BOLD_signal   : array of BOLD signal (x,y,t)
%                         size(BOLD_signal) = [length(y), length(x), length(t)]
%         deconvResponses    : array of deconvolved responses (x,y,t)
%                         size(deconvResponses.{}) = [length(y), length(x), length(t)]
%                         Possible fields are reconvBOLD, neural, neuroglial, 
%                         CBF, CBV, dHb, Wmode, Lmode, and Dmode.
%         x             : vector of distance along x
%         y             : vector of distance along y
%         t             : vector of time 
%         zoom_lim      : vector that limits the patch size of responses.
%                         zoom_lim(1,1) is the start location in x in mm scale
%                         zoom_lim(1,2) is the end location in x in mm scale
%                         zoom_lim(2,1) is the start location in y in mm scale
%                         zoom_lim(2,2) is the end location in y in mm scale
%         params        : instance of the class loadParameters of the 
%                         toolbox
%         normalization : 1 or 0
%                         Choose 1 if you want the responses normalized, 
%                         0 otherwise.
%         clim_factor   : dividing the maximum value of clim by this value
%         plot_what     : string saying which response(s) to plot.
%                         Possible inputs are 'BOLD', 'reconvBOLD', 'neural', 
%                         'neuroglial', 'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 
%                         'Dmode', 'all_w_BOLD', 'all_no_BOLD'.
%         movie_filename: filename of the saved movie
%         frame_rate    : frames per second of the movie
%         v1_boundary   : structure of pixel values of the boundary of v1.
%                         If this input is left empty, the code will still
%                         plot the responses but without the boundaries.
%                         Possible fields are 'x' and 'y' for the x and y
%                         pixel values of the boundary, respectively.
%
% Output: fig           : figure handle of the resulting movie
%
% James Pang, University of Sydney, 2016

%%
% find index of t=0
t0 = dsearchn(t', 0);

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

% Set up the movie.
myVideo = VideoWriter(movie_filename);       % Filename
myVideo.FrameRate = frame_rate;              % How many frames per second
myVideo.Quality = 100;   
open(myVideo); 

xmin_ind = dsearchn(x'*1e3, zoom_lim(1,1));
xmax_ind = dsearchn(x'*1e3, zoom_lim(1,2));
ymin_ind = dsearchn(y'*1e3, zoom_lim(2,1));
ymax_ind = dsearchn(y'*1e3, zoom_lim(2,2));

if strcmpi(plot_what, 'all_w_BOLD')
    width = 0.7; height = 0.098; initial_x = 0.19; initial_y = 0.87;
    y_factor = 1.1;

    fig = figure('Position', [200, 200, 180, 700], 'color','w');
    cmap = colormap_bluetored;
    
    for index=t0:length(t)
        for k=1:length(responses)-1
            if k==1
                data = real(BOLD_signal);
            else
                data = real(deconvResponses.(responses{k}));
            end
            
            clim = [min(min(min(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t))))), ...
                max(max(max(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t)))))]; 
            clim_max = max(abs(clim))/clim_factor;   

            subplot('Position', [initial_x initial_y-height*y_factor*(k-1) width height])
            
            imagesc(x(xmin_ind:xmax_ind)*1e3, y(ymin_ind:ymax_ind)*1e3, ...
                    data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,index))
            if nargin==13
                hold on;
                plot(x(v1_boundary.x)*1e3, y(v1_boundary.y)*1e3, 'k-', 'linewidth', 0.5)
                hold off;
            end    
            set(gca, 'fontSize', 11, 'xlim', [zoom_lim(1,1), zoom_lim(1,2)], ...
                'ylim', [zoom_lim(2,1), zoom_lim(2,2)], ...
                'xtick', [], 'ytick', [], 'YDir','normal');
            caxis([-clim_max, clim_max])
            colormap(cmap)
            colorbar
            
            if k==1
                title(sprintf('$t$ = %.2f s', t(index)), 'fontsize', 14, 'interpreter', 'latex')
            end
            ylabel(titles{k}, 'fontsize', 15)
        end
        
        frame = getframe(gcf); 
        writeVideo(myVideo, frame);
    end
elseif strcmpi(plot_what, 'all_no_BOLD')
    width = 0.7; height = 0.11; initial_x = 0.19; initial_y = 0.86;
    y_factor = 1.1;

    fig = figure('Position', [200, 200, 180, 600], 'color','w');
    cmap = colormap_bluetored;
    
    for index=t0:length(t)
        for k=2:length(responses)-1
            data = real(deconvResponses.(responses{k}));
            
            clim = [min(min(min(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t))))), ...
                max(max(max(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t)))))]; 
            clim_max = max(abs(clim))/clim_factor;   

            subplot('Position', [initial_x initial_y-height*y_factor*(k-2) width height])
            
            imagesc(x(xmin_ind:xmax_ind)*1e3, y(ymin_ind:ymax_ind)*1e3, ...
                    data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,index))
            if nargin==13
                hold on;
                plot(x(v1_boundary.x)*1e3, y(v1_boundary.y)*1e3, 'k-', 'linewidth', 0.5)
                hold off;
            end       
            set(gca, 'fontSize', 11, 'xlim', [zoom_lim(1,1), zoom_lim(1,2)], ...
                'ylim', [zoom_lim(2,1), zoom_lim(2,2)], ...
                'xtick', [], 'ytick', [], 'YDir','normal');
            caxis([-clim_max, clim_max])
            colormap(cmap)
            colorbar
            
            if k==2
                title(sprintf('$t$ = %.2f s', t(index)), 'fontsize', 14, 'interpreter', 'latex')
            end
            ylabel(titles{k}, 'fontsize', 15)
        end
        
        frame = getframe(gcf); 
        writeVideo(myVideo, frame);
    end
else
    width = 0.7; height = 0.65; initial_x = 0.19; initial_y = 0.15;

    fig = figure('Position', [200, 200, 180, 120], 'color','w');
    cmap = colormap_bluetored;
    
    response_index = find(strcmpi(responses, plot_what));
    if response_index==1
        data = real(BOLD_signal);
    else
        data = real(deconvResponses.(responses{response_index}));
    end
    clim = [min(min(min(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t))))), ...
                max(max(max(data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,t0:length(t)))))]; 
    clim_max = max(abs(clim))/clim_factor;   
    
    for index=t0:length(t)
        subplot('Position', [initial_x initial_y width height])

        imagesc(x(xmin_ind:xmax_ind)*1e3, y(ymin_ind:ymax_ind)*1e3, ...
                data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,index))
        if nargin==13
            hold on;
            plot(x(v1_boundary.x)*1e3, y(v1_boundary.y)*1e3, 'k-', 'linewidth', 0.5)
            hold off;
        end       
        set(gca, 'fontSize', 11, 'xlim', [zoom_lim(1,1), zoom_lim(1,2)], ...
            'ylim', [zoom_lim(2,1), zoom_lim(2,2)], ...
            'xtick', [], 'ytick', [], 'YDir','normal');
        caxis([-clim_max, clim_max])
        colormap(cmap)
        colorbar
        title(sprintf('$t$ = %.2f s', t(index)), 'fontsize', 14, 'interpreter', 'latex')
        ylabel(titles{response_index}, 'fontsize', 15)
        
        frame = getframe(gcf); 
        writeVideo(myVideo, frame);
    end
end

close(myVideo); % Saves the movie.
