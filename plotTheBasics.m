function plotTheBasics(mat_name, mat_path, selection, trim_start, trim_end, export)

% Basic Flight Data Visualization from MAT File
%
% DESCRIPTION: 
%   This script allows for the visualization of basic flight data
%   from a .mat file containing synchronized UAV log data. It supports two 
%   modes of operation:
%   1. Interactive Mode: Upon execution without any input arguments, the user 
%      is prompted to select a .mat file using a file explorer. The script 
%      then processes the selected data and visualizes it according to 
%      default selections.
%   2. Programmatic Mode: The script can be called with specific parameters 
%      to visualize data from a predefined .mat file without user interaction. 
%      This mode allows for automation and integration into larger analysis workflows.
%
% INPUTS: 
%   mat_name    - Name of the .mat file containing the flight data (optional)
%   mat_path    - Path to the directory containing the .mat file (optional)
%   selection   - Integer specifying the type of plot to visualize (optional)
%   trim_start  - Start time in seconds for data trimming (optional)
%   trim_end    - End time in seconds for data trimming (optional)
%   export      - True/False if the plot will be exported as a PDF 
%
% OUTPUTS:
%   Figures displaying the selected flight data visualizations
% 
% WRITTEN BY:
%   Garrett D. Asper
%   Virginia Tech
%   Email: garrettasper@vt.edu
%
% HISTORY:
%   03 JUL 2024 - Created and debugged, GDA
%
% NOMENCLATURE:                                     UNITS:
%   x         - inertial north position along i     (ft)
%   y         - inertial east position along j      (ft)    
%   z         - inertial north position along k     (ft)
%   u         - inertial frame velocity along i     (ft/s)
%   v         - inertial frame velocity along j     (ft/s)
%   w         - inertial frame velocity along k     (ft/s)
%   phi       - roll angle                          (deg)
%   theta     - pitch angle                         (deg)
%   psi       - yaw angle                           (deg)
%   p         - roll rate                           (deg/s)
%   q         - pitch rate                          (deg/s)
%   r         - yaw rate                            (deg/s)
%
% The following data are visualized based on the user's selection:
% (defined using the 'selection' variable)
%
%   1 - Plots 4 plots, includes 1. NED, 2. uvw, 3. phi theta psi, 4. pqr
%       
%   2 - phi theta r vs cmd
%   3 - Alt vs Alt_sp
%   4 - RPM vs RPM_cmd
%   5 - VV vs VV_sp
%
% Note: If running in interactive mode, after clicking "Run", a file explorer window 
%       will open, allowing the user to navigate to and select the desired .mat file for visualization.
%
% Ensure that all variables in the 'User Defined Data' section are defined
% before running the script

%% User Defined Data 

% ~~~~ Interactive Mode ~~~~ (Define these and then run the script)
if nargin == 0 || isempty(mat_name) || isempty(mat_path)
    clear; clc;
    % No input provided, Prompt the user to select the .mat file using the file explorer
    mat_path = strrep(mfilename('fullpath'), mfilename, '');
    [mat_name, mat_path] = uigetfile([mat_path, '*.mat'], 'MultiSelect', 'off');
    
    load(fullfile(mat_path, mat_name))

    % Default values if not provided
    selection = 4; % Example default value
    trim_start = 1; % Example default value
    trim_end = 0; % Example default value
    export = false;

% ~~~~ Programmatic Mode ~~~~ (Defined externally when the function is called)
elseif nargin < 6
    % In case not all parameters are provided, set default for missing ones
    if isempty(selection)
        error('No plot selection was made.')
    end
    if isempty(trim_end)
        trim_end = 0; % Default trim_end
    end
    if isempty(trim_start)
        trim_start = 0; % Default trim_start
    end
    load(fullfile(mat_path, mat_name))
elseif nargin == 6
    load(fullfile(mat_path, mat_name))
end

addpath('plotting_functions', 'utils'); % Add subfolders to the path

%% Format Setup
formatPlot() % Apply plot formatting

%% Time Vector Creation
time = sensor_combined.timestamp; % Import the synchronized time (duration array)
[time_seconds, start_idx, end_idx] = trimLogIdx(time, trim_start, trim_end);

%% Import Log Data
% Log data structure that contains the essentials for plotting 
log = createLogStructure(mat_name, mat_path);

% Trim log data based on specified time indices
log = trimLogData(log, start_idx, end_idx, time_seconds);

%% Call the Applicable Plotting Function

switch selection
    case 1
        exportTitle = '12 DOF Analysis';
        fig = stacked_12dof(seconds(time_seconds), log.all_12dof);
        set(get(gcf,'Children'),'FontSize', 12);
    case 2
        exportTitle = 'Innerloop Controller Performance';
        inner = [log.phi, log.theta, log.r];
        inner_sp = [log.phi_sp, log.theta_sp, log.r_sp];
        fig = innerLoop(time_seconds, inner, inner_sp);
        set(get(gcf,'Children'),'FontSize', 14);
    case 3
        exportTitle = 'Altitude Hold Performance';
        fig = alt_vs_sp(time_seconds, log.alt_ft, log.alt_ft_sp);
        set(get(gcf,'Children'),'FontSize', 14);
    case 4
        exportTitle = 'RPM Tracking';
        fig = RPM(time_seconds, log.rpm, log.rpm_sp);
        set(get(gcf,'Children'),'FontSize', 14);
    case 5
        exportTitle = 'Vertical Velocity vs Setpoint';
        fig = vv_vs_sp(time_seconds, log.vv_fps, log.vv_fps_sp);
        set(get(gcf,'Children'),'FontSize', 14);
    otherwise
        error('Invalid selection');
end

% Define the size of the figure when it is generated
set(gcf, 'Position', [0, 0, 1250, 600]); 

%% Export PDF
if export == true
    exportgraphics(fig,[exportTitle,'.pdf'],'ContentType','vector')
end

end