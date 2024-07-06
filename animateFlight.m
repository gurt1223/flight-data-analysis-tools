function animateFlight(ulog_name, ulog_path, freq, cutoffFrequency, trim_start, trim_end, mp4FileName)

% Basic Flight Animation from MAT File
%
% DESCRIPTION: 
%   This script generates an animation from ULOG data. It supports two 
%   modes of operation:
%   1. Interactive Mode: Without input arguments, it prompts the user to select
%      a ULOG file for conversion through a graphical interface and input 
%      parameters through dialog boxes.
%   2. Programmatic Mode: With provided parameters, it directly processes 
%      the specified ULOG file without user interaction.
%
% INPUTS: 
%   ulog_name       - Name of the ULog file (optional)
%   ulog_path       - Path to the directory containing the ULog file (optional)
%   freq            - Desired frequency for time synchronization (Hz, optional)
%   cutoffFrequency - Cutoff frequency for the low-pass filter (Hz, optional)
%   trim_start      - Start time in seconds for data trimming (optional)
%   trim_end        - End time in seconds for data trimming (optional)
%   mp4FileName     - Name of the output MP4 file (optional)
%
% OUTPUTS:
%   An MP4 video file visualizing the flight data
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

% Check if function is run without arguments
if nargin == 0
    clear; clc; clearvars -except videoObject; close all;
    % Prompt the user for parameters using an input dialog box
    prompt = {'Enter the frames per second of video:', ...
              'Enter the cutoff frequency for filtering (in Hz):', ...
              'Enter the start time (in seconds):', ...
              'Enter the end time (in seconds, enter 0 to see entire flight):', ...
              'Enter the name of the MP4 file (without extension):'};
    dlgtitle = 'Input Parameters';
    dims = [1 50];
    defaultInput = {'60', '10', '0', '0', 'myVideo'};
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

% Add the required files to the path
addpath('plotting_functions', 'utils','animate_flight'); % Add subfolders to the path

% Process the data from the mat file
processUlogTimeSynced(ulog_path, ulog_name, freq, cutoffFrequency)
mat_path = ulog_path;
mat_name = strrep(ulog_name, '.ulg', '.mat'); % Replace .ulg with .mat
load(fullfile(mat_path, mat_name))

% Create the log structure
log = createLogStructure(mat_name, mat_path);

% Time Vector Creation
time = sensor_combined.timestamp; % Import the synchronized time (note this is a duration array)
[time_seconds, start_idx, end_idx] = trimLogIdx(time, trim_start, trim_end);

% Trim the log data
log = trimLogData(log, start_idx, end_idx, time_seconds);

% Extract translation and rotation data
translations = [log.N, log.E, log.alt_ft];
rotations = log.quat;

% Aircraft options
scaleFactor = 3;
STL = 'quad.stl';

% Make the figure fairly large for good resolution
f = figure('Position', [50, 50, 1080, 1080], 'Visible', 'on');
f.Color = 'w';
box on;

% Check if the user pressed cancel or did not enter a filename
if isempty(mp4FileName)
    error('No filename entered. Exiting script.');
else
    fullMp4FileName = strcat(mp4FileName, '.mp4'); % Append the .mp4 extension
end
videoObject = VideoWriter(fullMp4FileName, 'MPEG-4');
videoObject.FrameRate = freq; % Match video frame rate with data sample rate

% Note: If you cancel the script mid-animation, you need to close the video object with close(videoObject)

% Create animation
tic;
animateObject(translations, rotations, scaleFactor, videoObject, STL, trim_start);
endTime = toc;

fprintf('Video rendering completed in %.2f seconds.\n', endTime);

end