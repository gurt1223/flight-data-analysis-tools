
% Description: This script allows the user to properly call the 
%   plotTheBasics function, which facilitates the visualization of basic 
%   flight data from a .mat file. By setting up input parameters such as the 
%   file name, path, selection type, and optional data trimming, users can 
%   easily generate visualizations either interactively or programmatically. 

clear; clc;

% Define the name of the .mat file containing the flight data
mat_name = '20240327_Flight2.mat';

% Define the path to the directory containing the .mat file
mat_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240327\';

% Specify the type of plot to visualize:
%   1 - NED, uvw, phi theta psi, pqr
%   2 - NED, uvw
%   3 - phi theta psi, pqr
%   4 - NED
%   5 - uvw
%   6 - phi theta psi
%   7 - pqr
%   8 - phi theta r vs cmd
%   9 - D vs D_cmd
%   10 - RPM vs RPM_cmd
%   11 - Plots 4 plots, includes 1. NED, 2. uvw, 3. phi theta psi, 4. pqr

selection = 8;

% Specify the start time in seconds for data trimming (optional)
trim_start = 1;

% Specify the end time in seconds for data trimming (optional)
trim_end = 0;

% Call the function to visualize the flight data based on the provided inputs
plotTheBasics(mat_name, mat_path, selection, trim_start, trim_end);
