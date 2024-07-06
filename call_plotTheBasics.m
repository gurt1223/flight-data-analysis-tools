% Example Call 'plotTheBasics.m' Function Script
%
% DESCRIPTION: 
%   This script allows the user to properly call the plotTheBasics function,
%   which facilitates the visualization of basic flight data from a .mat file.
%   By setting up input parameters such as the file name, path, selection type,
%   and optional data trimming, users can easily generate visualizations 
%   either interactively or programmatically.
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
%
%       1 - Plots 4 plots, includes 1. NED, 2. uvw, 3. phi theta psi, 4.
%       pqr
%       2 - phi theta r vs cmd
%       3 - Alt vs Alt_sp
%       4 - RPM vs RPM_cmd
%       5 - VV vs VV_sp

clear; clc;

% Define the name of the .mat file containing the flight data
mat_name = '20240327_Flight2.mat';

% Define the path to the directory containing the .mat file
mat_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240327\';

% Choose the selection number based on the desired plot
selection = 4;

% Specify the start time in seconds for data trimming (optional)
trim_start = 1;

% Specify the end time in seconds for data trimming (optional)
trim_end = 10;

% Select if the plot will be exported to a PDF
export = 'false'; % either 'true' or 'false'

% Call the function to visualize the flight data based on the provided inputs
plotTheBasics(mat_name, mat_path, selection, trim_start, trim_end, export);