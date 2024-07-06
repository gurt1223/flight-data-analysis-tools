% Example Call Process ULog File and Apply Time Synchronization
%
% DESCRIPTION:
%   This script processes a ULog file by applying time synchronization and
%   filtering. The ULog file is loaded from the specified path, and the 
%   data is synchronized at a specified frequency. A low-pass filter is 
%   applied to the data with a specified cutoff frequency.
%
% INPUTS:
%   ulog_name       - Name of the ULog file (string)
%   ulog_path       - Path to the directory containing the ULog file (string)
%   freq            - Desired frequency for time synchronization (Hz)
%   cutoffFrequency - Cutoff frequency for the low-pass filter (Hz)
%
% OUTPUTS:
%   Processed and synchronized ULog data saved to the specified directory
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

clear; clc;

% Define the name of the .ulg file containing the flight data
ulog_name = '20240327_Flight2.ulg';

% Define the path to the directory containing the .ulg file
ulog_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240327\';

% Define processing parameters
freq = 100; % Desired frequency for time synchronization [Hz]
cutoffFrequency = 10; % Cutoff frequency for the low-pass filter [Hz]

% Process ULog file with time synchronization and filtering
processUlogTimeSynced(ulog_path, ulog_name, freq, cutoffFrequency);