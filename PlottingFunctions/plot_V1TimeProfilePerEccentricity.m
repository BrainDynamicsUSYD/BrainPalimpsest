function fig = plot_V1TimeProfilePerEccentricity(responses_1D, eccentricity, ...
                                        eccentricity_interest, t, plot_what)
%% plot_V1TimeProfilePerEccentricity.m
%
% Plots the normalized time profiles of the responses in V1 for different
% eccentricities.
%
% Inputs: responses_1D          : array of average 1D responses in V1 BOLD (eccentricity, t)
%                                 size(responses_1D.{}) = [length(eccentricity), length(t)]
%                                 Possible fields are BOLD, reconvBOLD, 
%                                 neural, neuroglial, CBF, CBV, dHb, Wmode, 
%                                 Lmode, and Dmode.
%         eccentricity          : vector of eccentricity values
%         eccentricity_interest : vector of chosen eccentricity values that
%                                 we are interested to analyze
%         t                     : vector of time 
%         plot_what             : string saying which response to plot.
%                                 Possible inputs are 'BOLD', 'reconvBOLD', 
%                                 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', 
%                                 'Wmode', 'Lmode', 'Dmode', 'all_w_BOLD', 
%                                 'all_no_BOLD'.
%
% Output: fig                   : figure handle of the resulting plot
%
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%% 

titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode', 'reconv BOLD'};
responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode', 'reconvBOLD'};
colors = {'k', 'r', 'b', 'k', 'r', 'b', 'k', 'r', 'b'};
lines = {'-', '-', '-', '--', '--', '--', '-', '-', '-'};
widths = [2, 2, 2, 2, 2, 2, 1, 1, 1];  

eccentricity_interest = sort(eccentricity_interest);
eccentricity_interest_ind = dsearchn(eccentricity', eccentricity_interest');

if strcmpi(plot_what, 'all_w_BOLD')
    if length(eccentricity_interest) > 1
        fig = figure('Position', [200, 200, 400, 120*length(eccentricity_interest)]);
        
        for i = 1:length(eccentricity_interest)
            minimum = [];

            sub = subplot(length(eccentricity_interest), 1, i);
            pos = get(gca, 'Position');
            delete(sub)
            subplot('Position', pos+[0.01,0,-0.2,0])

            hold on;
            for j=1:length(responses)-1
                data = real(responses_1D.(responses{j}));
                data_plot = data(:, eccentricity_interest_ind(i))'/max(data(:, eccentricity_interest_ind(i)));
                minimum = min([minimum, min(data_plot)]);

                plot(t, data_plot, 'Color', colors{j}, 'LineStyle', lines{j}, 'LineWidth', widths(j))
            end
            hold off;
            if i==length(eccentricity_interest)
                xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            end
            set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [minimum, 1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(minimum):0.5:1)

            text(-0.17*(t(end)-t(1)), minimum+0.05*(1-minimum), ...
                ['$E = $', num2str(eccentricity_interest(i)), '$^\circ$'], ...
                'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
            if i==1
                leg = legend(titles{1:length(responses)-1});
                set(leg, 'FontSize', 12, 'Orientation','Vertical',...
                    'box','off','Position',get(gca, 'Position')+[0.68,-0.24/length(eccentricity_interest),-0.5,0])
            end
        end
    else
        fig = figure('Position', [200, 200, 400, 120*3]);
        
        minimum = [];
        for j=1:length(responses)-1
            data = real(responses_1D.(responses{j}));
            data_plot = data(:, eccentricity_interest_ind)'/max(data(:, eccentricity_interest_ind));
            minimum = min([minimum, min(data_plot)]);
        end
            
        for i=1:3
            sub = subplot(3, 1, i);
            pos = get(gca, 'Position');
            delete(sub)
            subplot('Position', pos+[0.01,0,-0.2,0])

            hold on;
            for j=((1:3)+3*(i-1))
                data = real(responses_1D.(responses{j}));
                data_plot = data(:, eccentricity_interest_ind)'/max(data(:, eccentricity_interest_ind));

                plot(t, data_plot, 'Color', colors{j}, 'LineStyle', lines{j}, 'LineWidth', widths(j))
            end
            hold off;
            if i==3
                xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            end
            set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [minimum, 1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(minimum):0.5:1)
            if i==2
                text(-0.17*(t(end)-t(1)), minimum+0.05*(1-minimum), ...
                    ['$E = $', num2str(eccentricity_interest), '$^\circ$'], ...
                    'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
            end
            
            leg = legend(titles{(1:3)+3*(i-1)});
            set(leg, 'FontSize', 12, 'Orientation','Vertical',...
                'box','off','Position',get(gca, 'Position')+[0.68,0.15,-0.5,-0.2])
        end
    end
    
elseif strcmpi(plot_what, 'all_no_BOLD')
    if length(eccentricity_interest) > 1
        fig = figure('Position', [200, 200, 400, 120*length(eccentricity_interest)]);
        
        for i = 1:length(eccentricity_interest)
            minimum = [];

            sub = subplot(length(eccentricity_interest), 1, i);
            pos = get(gca, 'Position');
            delete(sub)
            subplot('Position', pos+[0.01,0,-0.2,0])

            hold on;
            for j=2:length(responses)-1
                data = real(responses_1D.(responses{j}));
                data_plot = data(:, eccentricity_interest_ind(i))'/max(data(:, eccentricity_interest_ind(i)));
                minimum = min([minimum, min(data_plot)]);

                plot(t, data_plot, 'Color', colors{j}, 'LineStyle', lines{j}, 'LineWidth', widths(j))
            end
            hold off;
            if i==length(eccentricity_interest)
                xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            end
            set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [minimum, 1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(minimum):0.5:1)

            text(-0.17*(t(end)-t(1)), minimum+0.05*(1-minimum), ...
                ['$E = $', num2str(eccentricity_interest(i)), '$^\circ$'], ...
                'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
            if i==1
                leg = legend(titles{2:length(responses)-1});
                set(leg, 'FontSize', 12, 'Orientation','Vertical',...
                    'box','off','Position',get(gca, 'Position')+[0.68,-0.24/length(eccentricity_interest),-0.5,0])
            end
        end
    else
        fig = figure('Position', [200, 200, 400, 120*3]);
        minimum = [];
        for j=1:length(responses)-1
            data = real(responses_1D.(responses{j}));
            data_plot = data(:, eccentricity_interest_ind)'/max(data(:, eccentricity_interest_ind));
            minimum = min([minimum, min(data_plot)]);
        end
        
        for i=1:3
            sub = subplot(3, 1, i);
            pos = get(gca, 'Position');
            delete(sub)
            subplot('Position', pos+[0.01,0,-0.2,0])

            hold on;
            for j=((1:3)+3*(i-1))
                if j~=1
                    data = real(responses_1D.(responses{j}));
                    data_plot = data(:, eccentricity_interest_ind)'/max(data(:, eccentricity_interest_ind));

                    plot(t, data_plot, 'Color', colors{j}, 'LineStyle', lines{j}, 'LineWidth', widths(j))
                end
            end
            hold off;
            if i==3
                xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            end
            set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [minimum, 1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(minimum):0.5:1)
            if i==2
                text(-0.17*(t(end)-t(1)), minimum+0.05*(1-minimum), ...
                    ['$E = $', num2str(eccentricity_interest), '$^\circ$'], ...
                    'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
            end
            
            if i==1
                leg = legend(titles{2:3});
            else
                leg = legend(titles{(1:3)+3*(i-1)});
            end
            set(leg, 'FontSize', 12, 'Orientation','Vertical',...
                'box','off','Position',get(gca, 'Position')+[0.68,0.15,-0.5,-0.2])
        end
    end
    
else
    if length(eccentricity_interest) > 1
        fig = figure('Position', [200, 200, 350, 120*length(eccentricity_interest)]);

        for i = 1:length(eccentricity_interest)
            sub = subplot(length(eccentricity_interest), 1, i);
            pos = get(gca, 'Position');
            delete(sub)
            subplot('Position', pos+[0.05,0,-0.02,0])

            response_index = find(strcmpi(responses, plot_what));
            data = real(responses_1D.(responses{response_index}));
            data_plot = data(:, eccentricity_interest_ind(i))'/max(data(:, eccentricity_interest_ind(i)));

            plot(t, data_plot, 'ko-', 'LineWidth', 2, 'MarkerSize', 3)

            if i==length(eccentricity_interest)
                xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            end
            set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [min(data_plot)*1.1, 1.1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(min(data_plot)):0.5:1)

            text(-0.17*(t(end)-t(1)), min(data_plot)+0.05*(1-min(data_plot)), ...
                ['$E = $', num2str(eccentricity_interest(i)), '$^\circ$'], ...
                'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
            if i==1
                title(titles{response_index}, 'fontsize',15);
            end
        end
    else
        fig = figure('Position', [200, 200, 350, 120]);
        
        sub = subplot(1, 1, 1);
        pos = get(gca, 'Position');
        delete(sub)
        subplot('Position', pos+[0.03,0.2,-0.2,-0.2])
        
        response_index = find(strcmpi(responses, plot_what));
        data = real(responses_1D.(responses{response_index}));
        data_plot = data(:, eccentricity_interest_ind)'/max(data(:, eccentricity_interest_ind));

        plot(t, data_plot, 'k-', 'LineWidth', 2)
        xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
        set(gca, 'FontSize', 12, 'XLim', [t(1), t(end)], ...
                'YLim', [min(data_plot)*1.1, 1.1], 'XTick', t(1):5:t(end), ...
                'YTick', floor(min(data_plot)):0.5:1)

        text(-0.17*(t(end)-t(1)), min(data_plot)+0.05*(1-min(data_plot)), ...
            ['$E = $', num2str(eccentricity_interest), '$^\circ$'], ...
            'fontsize', 18, 'fontweight', 'b', 'rotation', 90, 'interpreter', 'latex')
        
        leg = legend(titles{response_index});
        set(leg, 'FontSize', 12, 'Orientation','Vertical',...
            'box','off','Position',get(gca, 'Position')+[0.68,0.15,-0.5,-0.2])
    end
end
    