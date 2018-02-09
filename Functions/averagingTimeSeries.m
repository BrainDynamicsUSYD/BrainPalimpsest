function [avg_grid] = averagingTimeSeries(scanNo, grid_data)
%% AveragingTimeSeries.m
%
% Calculates the time average of grid_data across different different 
% cycles of the scan. This reduces the time series into a single period 
% with length of cycles_length.
%
% Inputs: scanNo        : scan number to be processed
%                         Possible fields are 
%                         1 for the expanding ring data
%                         2, 3, ..., 11 for the expanding arc data.
%         grid_data     : matrix of data
%
% Output: avg_grid      : matrix of average of grid_data
% 
% Original: James Pang, University of Sydney, 2017
% Version 1.2: James Pang, University of Sydney, Jan 2018

%% 

% a list to sort out the timing files and relate them to the scan
timingFileList = [6, 7, 9, 10, 11, 12, 2, 3, 4, 5];

% finding the size of grid_data data
[nx, ny, ~] = size(grid_data);    
    
if scanNo == 1
    cycles = 11;
    cycles_length = 15;
    
    relevant_index = reshape(16:180, [cycles_length, cycles])';
else
    load(['Data/ExpandingRingAndExpandingArc/TimingFiles/Scan' num2str(timingFileList(scanNo - 1)) '.mat'], ...
        'start_trial_pulses');
    
    cycles = length(start_trial_pulses);
    cycles_length = start_trial_pulses(2) - start_trial_pulses(1) - 1;
    
    relevant_index = zeros(cycles, cycles_length);
    for k=1:length(start_trial_pulses)
        relevant_index(k, :) = start_trial_pulses(k):cycles_length+start_trial_pulses(k)-1;
    end      
end

relevant_index = sort(relevant_index(:));  

relevant_data = grid_data(:, :, relevant_index);
nt = size(relevant_data, 3);

% reshaping and doing the average
reshaped_relevant_data = reshape(relevant_data, [nx*ny, nt]);
newData = reshape(reshaped_relevant_data, [nx*ny, cycles_length, cycles]);
newData = mean(newData, 3);

% reshaping and assigning the averaged data to a new variable
avg_grid = reshape(newData, [nx, ny, cycles_length]);
