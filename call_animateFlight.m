% Example Call 'animateFlight.m' Function Script
%
% DESCRIPTION: 
%   This script allows the user to properly call the animateFlight function,
%   which generates an animation from ULOG data. By setting up input parameters
%   such as the file name, path, frequency, cutoff frequency, start and end times,
%   and output file name, users can easily create a flight animation either 
%   interactively or programmatically.
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

% Define the name of the .ulg file containing the flight data
ulog_name = '19_40_36.ulg';

% Define the path to the directory containing the .ulg file
ulog_path = 'C:\Users\gasper\OneDrive - NASA\Documents\IMPACT Vehicle\Aircraft Logbook\20240702\';

freq = 60; % [Hz] (also results in frames per second of video)
cutoffFrequency = 10; % [Hz] for filtering angular rates
trim_start = 0; % Start time of the video corresponding to flight log time (seconds)
trim_end = 10; % End time of the video corresponding to flight log time (seconds)
               % Leave 0 for full flight log
mp4FileName = 'Test1'; % Name of the visualization video file (do not append with .mp4)

% Call the function
animateFlight(ulog_name, ulog_path, freq, cutoffFrequency, trim_start, trim_end, mp4FileName)
