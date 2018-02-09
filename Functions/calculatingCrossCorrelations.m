function calculatingCrossCorrelations(hemisphere, resolution, what_correlation)
%% calculatingCrossCorrelations.m
%
% Calculates the cross correlations of the quantities.
%
% Inputs: hemisphere        : string of hemisphere
%                             Possible fields are lh for left hemisphere and 
%                             rh for right hemisphere.
%         resolution        : x and y spatial resolution in mm
%         what_correlation  : string of quantity to which other quantities
%                             will be correlated with.
%                             Possible inputs are 'visualStimulus', 'BOLD', 
%                             'reconvBOLD', 'neural', 'neuroglial', 'CBF', 
%                             'CBV', 'dHb', 'Wmode', 'Lmode', and 'Dmode'.
%
% Output: none
% 
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%
% hemisphere = 'lh';
% resolution = 0.2;
% what_correlation = 'BOLD';

dt = 2;
min_lag_resolution = 5;
shifts = min_lag_resolution:-1:-min_lag_resolution;

t_lags = -shifts*dt;
t_lags_interp = min(t_lags):0.1:max(t_lags);

%% Mean cross correlations

quantities = {'visualStimulus', 'BOLD', 'reconvBOLD', 'neural', 'neuroglial', ...
             'CBF', 'CBV', 'dHb', 'Wmode', 'Lmode', 'Dmode'};

% Calculating cross correlations for all scans and responses
for scanNo = 1:11
    filename1 = ['Data/ExpandingRingAndExpandingArc/GriddedMatFiles/', ...
                  hemisphere, '.Scan', num2str(scanNo), '_resolution=', ...
                  num2str(resolution), '.mat'];
    filename2 = ['Data/ExpandingRingAndExpandingArc/VisualStimulus/', ...
                  hemisphere, '.Scan', num2str(scanNo), '_VisualStimulus_resolution=', ...
                  num2str(resolution), '.mat'];
    if scanNo == 1
        load(filename1, 'reordered_avg_BOLD_signal', 'reordered_deconvResponses_avg');
        avg_BOLD_signal = reordered_avg_BOLD_signal;
        deconvResponses_avg = reordered_deconvResponses_avg;
    else
        load(filename1, 'avg_BOLD_signal', 'deconvResponses_avg');
    end
    
    load(filename2, 'visualStimulus_smooth');

    what_correlation_index = find(strcmpi(quantities, what_correlation));
    if what_correlation_index == 1
        data1 = visualStimulus_smooth;
    elseif what_correlation_index == 2
        data1 = avg_BOLD_signal;
    else
        data1 = deconvResponses_avg.(quantities{what_correlation_index});
    end

    data1_diff = bsxfun(@minus, data1, mean(data1, 3));

    basis = avg_BOLD_signal(:,:,1)./avg_BOLD_signal(:,:,1);

    for i = 1:length(quantities)
        if i == 1
            data2 = visualStimulus_smooth;
        elseif i == 2
            data2 = avg_BOLD_signal;
        else
            data2 = deconvResponses_avg.(quantities{i});
        end

        cross_correlation = zeros(1, length(shifts));
        for lag = 1:length(shifts)
            data2_shifted = data2(:,:,circshift((1:size(avg_BOLD_signal, 3))', shifts(lag))');

            data2_diff = bsxfun(@minus, data2_shifted, mean(data2_shifted, 3));

            r = sum(data1_diff.*data2_diff, 3)./sqrt(sum(data1_diff.^2, 3).*sum(data2_diff.^2, 3));
            result = basis;
            result(~isnan(basis)) = r(~isnan(basis));

            cross_correlation(lag) = mean(result(~isnan(result)));
        end

        correlations.(quantities{i})(scanNo, :) = cross_correlation; 
    end
end

%% Interpolating actual cross correlations
for i = 1:length(quantities)
    for scanNo = 1:3%11
        interp_correlations.(quantities{i})(scanNo, :) = interp1(t_lags, ...
                correlations.(quantities{i})(scanNo, :), t_lags_interp, 'spline');
    end 
end

%% Taking the mean and std across all scans
for i = 1:length(quantities)
    mean_correlations.(quantities{i}) = mean(correlations.(quantities{i}),1);
    std_correlations.(quantities{i}) = std(correlations.(quantities{i}), 0, 1);
end

%% Interpolating mean cross correlations
for i = 1:length(quantities)
    interp_mean_correlations.(quantities{i}) = interp1(t_lags, ...
                mean_correlations.(quantities{i}), t_lags_interp, 'spline');
end

%%             
filename3 = ['Data/ExpandingRingAndExpandingArc/CorrelationMatFiles/',...
             hemisphere,'.CrossCorrelations_',what_correlation,'.mat'];
save(filename3, 'correlations', 'interp_correlations', 'mean_correlations', ...
                'std_correlations', 'interp_mean_correlations', 't_lags', ...
                't_lags_interp')
