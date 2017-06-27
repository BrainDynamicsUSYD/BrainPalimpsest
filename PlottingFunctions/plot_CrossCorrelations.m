function fig = plot_CrossCorrelations(hemisphere, what_correlation)
%% plot_CrossCorrelations.m
%
% Plots the average cross-correlation of the quantity corresponding to the
% input what_correlation and the other quantities.
%
% Inputs: hemisphere        : string of hemisphere
%                             Possible fields are lh for left hemisphere 
%                             and rh for right hemisphere.
%         what_correlation  : string of quantity to which other
%                             responses will be correlated with.
%                             Possible inputs are 'visualStimulus', 'BOLD', 
%                             'reconvBOLD', 'neural', 'neuroglial', 'CBF', 
%                             'CBV', 'dHb', 'Wmode', 'Lmode', and 'Dmode'.
%
% Output: fig               : figure handle of the resulting plot
% 
% James Pang, University of Sydney, 2017

%% Loading the correlation file

filename = ['Data/ExpandingRingAndExpandingArc/CorrelationMatFiles/' ,...
            hemisphere,'.CrossCorrelations_',what_correlation,'.mat'];
load(filename, 'correlations', 'interp_correlations', 'mean_correlations', ...
               'std_correlations', 'interp_mean_correlations', 't_lags', ...
               't_lags_interp')

%% Plotting the results

responses = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
             'Wmode', 'Lmode', 'Dmode'};
titles = {'BOLD', 'neural', 'neuroglial', 'CBF', 'CBV', 'dHb', ...
          '{\it W} mode', '{\it L} mode', '{\it D} mode'};
      
what_correlation_index = find(strcmpi(responses, what_correlation));

plot_order = [what_correlation_index, 1:(what_correlation_index-1), ...
              (what_correlation_index+1):length(responses)];

width = 0.7; height = 0.1; initial_x = 0.23; initial_y = 0.88;
y_factor = 1;

dt = mean(diff(t_lags));

fig = figure('Position', [200, 200, 250, 800]);

for j=1:length(responses)
    subplot(length(responses), 1, j, 'Parent', fig, 'Position', [initial_x initial_y-height*y_factor*(j-1) width height])

    data_mean = mean_correlations.(responses{plot_order(j)});
    data_std = std_correlations.(responses{plot_order(j)});
    data_interp = interp_mean_correlations.(responses{plot_order(j)});
    [~, peak_ind] = max(abs(data_interp));

    errorbar(t_lags, data_mean, data_std, 'k.', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
    hold on;
    plot(t_lags_interp, data_interp, 'k-')
    plot([t_lags_interp(peak_ind), t_lags_interp(peak_ind)], [-1, 1], 'r--', 'LineWidth', 1.5) 
    hold off;
    if data_interp(peak_ind) > 0
        text(t_lags_interp(peak_ind)+0.5, -0.5, num2str(t_lags_interp(peak_ind)), 'FontSize', 10, 'FontWeight', 'b');
    else
        text(t_lags_interp(peak_ind)+0.5, 0.5, num2str(t_lags_interp(peak_ind)), 'FontSize', 10, 'FontWeight', 'b');
    end
    text(-10, 0.75, titles{plot_order(j)}, 'FontSize', 12, 'FontWeight', 'b')
    set(gca, 'FontSize', 12, 'XLim', [t_lags(1)*1.05, t_lags(end)*1.05], ...
            'YLim', [-1, 1], 'XTick', t_lags(1):dt:t_lags(end), ...
            'YTick', -1:0.5:1, 'YTickLabel', {'','-0.5','0','0.5',''}, ...
            'ticklength', get(gca,'ticklength'));%'TickLength', [0.003, 0.01])
    if j==9
        xlabel('lag (s)', 'fontsize', 15, 'interpreter', 'latex')
    elseif j==5
        ylabel(['cross-correlation with ', titles{what_correlation_index}], 'fontsize', 15, 'interpreter', 'latex')
        set(gca, 'XTickLabel', {});
    else
        set(gca, 'XTickLabel', {});
    end

end
  