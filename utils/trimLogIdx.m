function [time_seconds, start_idx, end_idx] = trimLogIdx(time, trim_start, trim_end)

% Trim Log Indices Based on Start and End Times
%
% DESCRIPTION:
%   This function calculates the time vector in seconds relative to the 
%   start time and determines the start and end indices for trimming based 
%   on specified start and end times.
%
% INPUTS:
%   time        - Time vector as a duration array
%   trim_start  - Start time in seconds for trimming
%   trim_end    - End time in seconds for trimming
%
% OUTPUTS:
%   time_seconds - Trimmed time vector in seconds
%   start_idx    - Start index for trimming
%   end_idx      - End index for trimming
% 
% WRITTEN BY:
%   Garrett D. Asper
%   Virginia Tech
%   Email: garrettasper@vt.edu
%
% HISTORY:
%   03 JUL 2024 - Created and debugged, GDA
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND

    % Calculate the time vector in seconds relative to the start time
    time_seconds = seconds(time - time(1));
    
    % Determine the end time based on the trim_end parameter
    if trim_end == 0
        end_time = time_seconds(end);
    else
        % Determine the range of times to be plotted
        end_time = time_seconds(end) - trim_end;
    end
    
    % Find the indices corresponding to the start and end times
    start_idx = find(time_seconds >= trim_start, 1);
    end_idx = find(time_seconds <= end_time, 1, 'last');
    time_seconds = time_seconds(start_idx:end_idx);
end
