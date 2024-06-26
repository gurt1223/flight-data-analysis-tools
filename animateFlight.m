function animateFlight(ulog_name, ulog_path,freq, cutoffFrequency, trim_start, trim_end, mp4FileName)
%generateAnimation Generates an animation from ULOG data.
% This function can be run in an interactive manner or by passing parameters directly.

% Check if function is run without arguments
if nargin == 0
    clear;clc; clearvars -except videoObject; close all;
    % Prompt the user for parameters using an input dialog box
    prompt = {'Enter the frames per second of video:', ...
              'Enter the cutoff frequency for filtering (in Hz):', ...
              'Enter the start time (in seconds):', ...
              'Enter the end time (in seconds, enter 0 to see entire flight):', ...
              'Enter the name of the MP4 file (without extension):'};
    dlgtitle = 'Input Parameters';
    dims = [1 50];
    defaultInput = {'60', '10','0', '0', 'myVideo'};
    userInput = inputdlg(prompt, dlgtitle, dims, defaultInput);
    
    % Check if the user pressed cancel
    if isempty(userInput)
        error('No input provided. Exiting script.');
    end
    
    % Parse user input
    freq = str2double(userInput{1});
    cutoffFrequency = str2double(userInput{2});
    trim_start = str2double(userInput{3});
    trim_end = str2double(userInput{4});
    mp4FileName = userInput{5};

    % Interactively select a ULOG file and process it
    ulog_path = strrep(mfilename('fullpath'), mfilename, '');
    [ulog_name, ulog_path] = uigetfile([ulog_path, '*.ulg'], 'MultiSelect', 'off');
end

% Process ulg into mat
% processUlogTimeSynced(ulog_path, ulog_name, freq, cutoffFrequency);

% Process the data from the mat file
mat_path = ulog_path;
mat_name = strrep(ulog_name, '.ulg', '.mat'); % Replace .ulg with .mat
load(fullfile(mat_path, mat_name))

% Time Vector Creation
time = sensor_combined.timestamp; % Import the synchronized time (note this is a duration array)
time_seconds = seconds(time - time(1)); % Create time vector in seconds

if trim_end == 0
    end_time = time_seconds(end);
else
    % Determine the range of times to be plotted
    %end_time_ = time_seconds(end) - trim_end;
    end_time = time_seconds(end) - (time_seconds(end) - trim_end);
end

% Find the indices corresponding to the start and end times
start_index = find(time_seconds >= trim_start, 1);
end_index = find(time_seconds <= end_time, 1, 'last');

ts = seconds(trim_start);

% Import relevant translation and rotation data
N = vehicle_local_position.x;
E = vehicle_local_position.y;
D = vehicle_local_position.z;
q = vehicle_attitude.q;
q1 = q(:,1); q2 = q(:,2); q3 = q(:,3); q4 = q(:,4);

% trim
N = N(start_index:end_index);
E = E(start_index:end_index);
D = D(start_index:end_index);
q1 = q1(start_index:end_index);
q2 = q2(start_index:end_index);
q3 = q3(start_index:end_index);
q4 = q4(start_index:end_index);


translations = [N, E, D]; % Example, replace with actual variable if different
rotations = [q1, q2, q3, q4]; % Convert quaternions to Euler angles if needed, replace with actual variables

% aircraft options
scaleFactor = 3;
STL = 'quad.stl';

% Make the figure fairly large for good resolution
% Make the figure fairly large for good resolution and invisible
f = figure('Position',[50 50 1080 1080], 'Visible', 'on');
f.Color = 'w';
box on

% Check if the user pressed cancel or did not enter a filename
if isempty(mp4FileName)
    error('No filename entered. Exiting script.');
else
    fullMp4FileName = strcat(mp4FileName, '.mp4'); % Append the .mp4 extension
end
videoObject = VideoWriter(fullMp4FileName,'MPEG-4');
videoObject.FrameRate = freq; % Match video frame rate with data sample rate
% Note: If you cancel the script mid-animation, you need to close the video
% object with close(videoObject)

% Create animation
tic
animateObject(translations, rotations, scaleFactor, videoObject, STL, ts);
endTime = toc;

fprintf('Video rendering completed in %.2f seconds.\n', endTime);

end