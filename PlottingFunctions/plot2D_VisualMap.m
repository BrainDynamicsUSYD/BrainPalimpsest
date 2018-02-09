function fig = plot2D_VisualMap(visualStimulus, x, y, t, ...
                    tslice, zoom_lim, v1_boundary)
%% plot2D_VisualMap.m
%
% Plots the 2D retinotopic map of the visual stimulus at different time slices.
%
% Inputs: visualStimulus: array of visual stimulus (x,y,t)
%                         size(visualStimulus) = [length(y), length(x), length(t)]
%         x             : vector of distance along x
%         y             : vector of distance along y
%         t             : vector of time 
%         tslice        : vector of time slices. Make sure that
%                         length(tslice)<=5 for now to fit in the
%                         predefined figure size.
%         zoom_lim      : vector that limits the patch size of responses.
%                         zoom_lim(1,1) is the start location in x in mm scale
%                         zoom_lim(1,2) is the end location in x in mm scale
%                         zoom_lim(2,1) is the start location in y in mm scale
%                         zoom_lim(2,2) is the end location in y in mm scale
%         v1_boundary   : structure of pixel values of the boundary of v1.
%                         If this input is left empty, the code will still
%                         plot the responses but without the boundaries.
%                         Possible fields are 'x' and 'y' for the x and y
%                         pixel values of the boundary, respectively.
%
% Output: fig           : figure handle of the resulting plot
%
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%
% find index of t=0
t0 = dsearchn(t', 0);

% time snapshot indices
tslice_ind = dsearchn(t', tslice');

xmin_ind = dsearchn(x'*1e3, zoom_lim(1,1));
xmax_ind = dsearchn(x'*1e3, zoom_lim(1,2));
ymin_ind = dsearchn(y'*1e3, zoom_lim(2,1));
ymax_ind = dsearchn(y'*1e3, zoom_lim(2,2));

width = 0.18; height = 0.67; initial_x = 0.01; initial_y = 0.12;
x_factor = 1.1; y_factor = 1.1;

fig = figure('Position', [200, 200, 450, 94]);
cmap = colormap_bluetored;

data = visualStimulus;

clim = [min(min(min(data(:,:,tslice_ind)))), max(max(max(data(:,:,tslice_ind))))];
clim_max = max(abs(clim));
for j=1:length(tslice)
    subplot('Position', [initial_x+width*x_factor*(j-1) initial_y-height*y_factor*0 width height])
    imagesc(x(xmin_ind:xmax_ind)*1e3, y(ymin_ind:ymax_ind)*1e3, ...
                    data(ymin_ind:ymax_ind,xmin_ind:xmax_ind,tslice_ind(j)))
    if nargin==7
        hold on;
        plot(x(v1_boundary.x)*1e3, y(v1_boundary.y)*1e3, 'k-', 'linewidth', 0.5)
        hold off;
    end
    set(gca, 'fontSize', 11, 'xlim', [zoom_lim(1,1), zoom_lim(1,2)], ...
            'ylim', [zoom_lim(2,1), zoom_lim(2,2)], ...
            'xtick', [], 'ytick', [], 'YDir','normal');
    title(['$t =$ ',num2str(tslice(j)), ' s'], 'fontsize', 14, 'interpreter', 'latex')
    caxis([-clim_max, clim_max])
    colormap(cmap)
end