
function plotTheBasics(mat_name, mat_path, selection, trim_start, trim_end)

% Basic Flight Data Visualization from MAT File
%
% (C) 2024 Garrett D. Asper <garrettasper@vt.edu>
% Analysis and Visualization of Flight Data
% 
% Description: This script allows for the visualization of basic flight data
% from a .mat file containing synchronized UAV log data. It supports two modes of operation:
%   1. Interactive Mode: Upon execution without any input arguments, the user 
%      is prompted to select a .mat file using a file explorer. The script 
%      then processes the selected data and visualizes it according to 
%      default selections.
%   2. Programmatic Mode: The script can be called with specific parameters 
%      to visualize data from a predefined .mat file without user interaction. 
%      This mode allows for automation and integration into larger analysis workflows.
%
% MATLAB Requirements: Matlab R2023a or newer
%
% Nomenclature:
%       x         - inertial north position along i
%       y         - inertial east position along j
%       z         - inertial north position along k 
%       u         - inertial frame velocity along i
%       v         - inertial frame velocity along j
%       w         - inertial frame velocity along k
%       phi       - roll angle
%       theta     - pitch angle
%       psi       - yaw angle
%       p         - roll rate
%       q         - pitch rate
%       r         - yaw rate
%
% The following data are visualized based on the user's selection:
% (defined using the 'selection' variable)
%
%       1 - NED, uvw, phi theta psi, pqr
%       2 - NED, uvw
%       3 - phi theta psi, pqr
%       4 - NED
%       5 - uvw
%       6 - phi theta psi
%       7 - pqr
%       8 - phi theta r vs cmd
%       9 - D vs D_cmd
%       10 - RPM vs RPM_cmd
%       11 - Plots 4 plots, includes 1. NED, 2. uvw, 3. phi theta psi, 4.
%       pqr
%
% Input(s): 
%       mat_name (optional): Name of the .mat file containing the flight data
%       mat_path (optional): Path to the directory containing the .mat file
%       selection (optional): Integer specifying the type of plot to visualize
%       trim_start (optional): Start time in seconds for data trimming
%       trim_end (optional): End time in seconds for data trimming
%
% Output(s): 
%       Figures displaying the selected flight data visualizations
%
% Usage:
%       Interactive Mode: plotTheBasics
%       Programmatic Mode: plotTheBasics('flightData.mat', 'C:/data', 8, 0, 180)
%
% Note: If running in interactive mode, after clicking "Run", a file explorer window 
%       will open, allowing the user to navigate to and select the desired .mat file for visualization.
%
% Ensure that all variables in the 'User Defined Data' section are defined
% before running the script
%
% Exporting:
%       To export as a pdf, run the following:
%       exportgraphics(fig,'fig_1.pdf','ContentType','vector')

%% User Defined Data 

% ~~~~ Interactive Mode ~~~~ (Define these and then run the script)
if nargin == 0 || isempty(mat_name) || isempty(mat_path)
    clear; clc;
    % No input provided, Prompt the user to select the .mat file using the file explorer
    mat_path = strrep(mfilename('fullpath'), mfilename, '');
    [mat_name, mat_path] = uigetfile([mat_path, '*.mat'], 'MultiSelect', 'off');
    
    load(fullfile(mat_path, mat_name))

    % Default values if not provided
    selection = 10; % Example default value
    trim_start = 1; % Example default value
    trim_end = 0; % Example default value
    export = true;

% ~~~~ 
elseif nargin < 5
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
elseif nargin == 5
    load(fullfile(mat_path, mat_name))
end

%% Format Setup
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',14);
set(0,'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on')   % Turn grid on for all plots
set(0, 'DefaultAxesYMinorTick', 'on');                      % Turn on minor tick marks for y-axis
set(0, 'DefaultAxesYTickMode', 'auto');                     % or 'manual' if you want to specify ticks manually
set(0,'DefaultAxesFontSize', 15)                            % Change the default figure axes font size
set(0, 'DefaultLineLineWidth', 1.1)                         % Set default plot line width
set(0, 'DefaultAxesLineWidth', 1.1)                         % Set default axes line width

% Define VT Colors
burnt_orange = [205/255 68/255 32/255];
chicago_maroon = [0.3882 0 0.1216]; 
sustainable_teal = [80/255 133/255 144/255];
% Define a new color order as an M-by-3 matrix of RGB values
newColorOrder = [burnt_orange;  % MATLAB default blue
                 chicago_maroon;  % MATLAB default orange
                 sustainable_teal;  % MATLAB default yellow
                 0.4940 0.1840 0.5560]; % MATLAB default purple

% Set the new color order globally for all plots
set(groot, 'DefaultAxesColorOrder', newColorOrder);

%% Time Vector Creation
time = sensor_combined.timestamp; % Import the synchronized time (note this is a duration array)
time_seconds = seconds(time - time(1)); % Create time vector in seconds

if trim_end == 0
    end_time = time_seconds(end);
else
    % Determine the range of times to be plotted
    end_timee = time_seconds(end) - trim_end;
    end_time = time_seconds(end) - end_timee;
end

% Find the indices corresponding to the start and end times
start_index = find(time_seconds >= trim_start, 1);
end_index = find(time_seconds <= end_time, 1, 'last');

%% Import Data from Retimed Log Topics
N = vehicle_local_position.x;
E = vehicle_local_position.y;
D = vehicle_local_position.z;

u = vehicle_local_position.vx;
v = vehicle_local_position.vy;
w = vehicle_local_position.vz;

q = vehicle_attitude.q;

eulXYZ = quat2eul(q,"ZYX");

% Trim the indices of the data
time_seconds = time_seconds(start_index:end_index);
phi  = eulXYZ(:,3)*(180/pi);
theta = eulXYZ(:,2)*(180/pi);
psi   = eulXYZ(:,1)*(180/pi);
p = vehicle_angular_velocity.xyz(:,1)*(180/pi);
q = vehicle_angular_velocity.xyz(:,2)*(180/pi);
r = vehicle_angular_velocity.xyz(:,3)*(180/pi);
p = p(start_index:end_index);
q = q(start_index:end_index);
r = r(start_index:end_index);
N = (N(start_index:end_index))*3.28084; % to feet
E = (E(start_index:end_index))*3.28084;
D = (D(start_index:end_index))*3.28084;
u = (u(start_index:end_index))*3.28084;
v = (v(start_index:end_index))*3.28084;
w = (w(start_index:end_index))*3.28084;
phi = phi(start_index:end_index);
theta = theta(start_index:end_index);
psi = psi(start_index:end_index);

% Check if 'custom_internal' is a timetable and if the new naming convention exists
if exist('custom_internal','var') && istimetable(custom_internal) && ismember('ctrl_thrust', custom_internal.Properties.VariableNames) 
    % New naming convention exists
    signals.ctrl_thrust = custom_internal.ctrl_thrust;
    signals.att_thrustcmds = custom_internal.att_thrustcmds;
    signals.z_sp = custom_internal.z_sp;
    signals.param_kp = custom_internal.param_kp;
    signals.param_ki = custom_internal.param_ki;
    signals.param_kd = custom_internal.param_kd;
    signals.param_kff = custom_internal.param_kff;
elseif exist('custom_internal','var')
    % Old naming convention fallback
    signals.ctrl_thrust(:,1) = custom_internal.signal_1;
    signals.ctrl_thrust(:,2) = custom_internal.signal_2;
    signals.ctrl_thrust(:,3) = custom_internal.signal_3;
    signals.ctrl_thrust(:,4) = custom_internal.signal_4;
    signals.att_thrustcmds(:,1) = custom_internal.signal_5;
    signals.att_thrustcmds(:,2) = custom_internal.signal_6;
    signals.att_thrustcmds(:,3) = custom_internal.signal_7;
    signals.att_thrustcmds(:,4) = custom_internal.signal_8;
    signals.z_sp = custom_internal.signal_9;
    signals.param_kp(:,1) = custom_internal.signal_10;
    signals.param_kp(:,2) = custom_internal.signal_11;
    signals.param_kp(:,3) = custom_internal.signal_12;
    signals.param_ki(:,1) = custom_internal.signal_13;
    signals.param_ki(:,2) = custom_internal.signal_14;
    signals.param_ki(:,3) = custom_internal.signal_15;
    signals.param_kd(:,1) = custom_internal.signal_16;
    signals.param_kd(:,2) = custom_internal.signal_17;
    signals.param_kd(:,3) = custom_internal.signal_18;
    signals.param_kff(:,1) = custom_internal.signal_19;
    signals.param_kff(:,2) = custom_internal.signal_20;
    signals.param_kff(:,3) = custom_internal.signal_21;
end

if selection ~= 10
    % Import and convert RC input
    ch1 = input_rc.values(:,1);
    ch1_trim = ch1(1)-1500;
    ch2 = input_rc.values(:,2);
    ch2_trim = ch2(10)-1500;
    ch4 = input_rc.values(:,4);
    ch4_trim = ch4(1)-1500;
    ch5 = input_rc.values(:,5);
    ch5 = ch5(start_index:end_index);
    ch4 = ch4(start_index:end_index);
    ch3 = input_rc.values(:,3);
    ch3 = ch3(start_index:end_index);
    ch3 = movmean(ch3,50);
    ch7 = input_rc.values(:,7);
    ch7 = ch7(start_index:end_index);
    
    % Assign an integer value to the flight mode switch channel
    condition1 = ch5 <= 1300;
    condition2 = ch5 <= 1600 & ch5 > 1300;
    condition3 = ch5 <= 2000 & ch5 > 1600;
    
    ch5_val = zeros(size(ch5));
    
    ch5_val(condition1) = 0; % altitude hold mode
    ch5_val(condition2) = 1; % stabilized mode
    ch5_val(condition3) = 2; % position mode

    % Do the same for Ch7 (PTI Activation)
    condition1 = ch7 <= 1300;
    condition2 = ch7 >= 1600;

    ch7_val(condition1) = 0; % PTI off
    ch7_val(condition2) = 1; % PTI on
    PTI_switch = find(abs(diff(ch7_val)) > 0); % catches when it's turned on

    % Stabilized mode commands
    roll_max_cmd = 60*(pi/180); % maximum cmd [rad]
    roll_trim = -(((ch1_trim))/400)*roll_max_cmd; % roll cmd [rad]    
    pitch_max_cmd = 60*(pi/180); % maximum cmd [rad]
    pitch_trim = (((ch2_trim))/400)*pitch_max_cmd; % roll cmd [rad]    
    yaw_max_cmd = -60*(pi/180); % maximum cmd [rad]
    yaw_trim = -(((ch4_trim))/400)*yaw_max_cmd; % yaw cmd [rad]

% minValue = 1100;
% maxValue = 1900;

    if exist('custom_internal', 'var')
        rollcmd_deg = movmean((signals.att_thrustcmds(:,1)*(180/pi))-(roll_trim*(180/pi)),10);
        pitchcmd_deg = movmean((signals.att_thrustcmds(:,2)*(180/pi))-(pitch_trim*(180/pi)),10);%(pitch_trim*(180/pi))),10);
        yawcmd_dps = movmean((signals.att_thrustcmds(:,3)*(180/pi)),10);
        throttlecmd = movmean((signals.att_thrustcmds(:,4)),5);
    
        % Trim (if requested)
        rollcmd_deg = rollcmd_deg(start_index:end_index);
        pitchcmd_deg = pitchcmd_deg(start_index:end_index);
        yawcmd_dps = yawcmd_dps(start_index:end_index);
        throttlecmd = throttlecmd(start_index:end_index);
        
        % Find max and min to set ylim
        rollcmd_max = max(max(rollcmd_deg,phi));
        rollcmd_min = min(min(rollcmd_deg,phi));
        pitchcmd_max = max(max(pitchcmd_deg,theta));
        pitchcmd_min = min(min(pitchcmd_deg,theta));
        yawcmd_max = max(max(yawcmd_dps,r));
        yawcmd_min = min(min(yawcmd_dps,r));
    end
end

if selection ~= 10
    if exist('custom_internal', 'var')
        z_sp_real = ((signals.z_sp)*3.28084);
        z_sp_real = z_sp_real(start_index:end_index);
    end
end

%% Determine which actuator_outputs name is used (compatibility btwn v1.13 and v1.14)
% After loading the .mat file
variableNames = who; % Get a list of variables in the workspace

% Check for the existence of actuator_outputs or actuator_outputs_1
if ismember('actuator_outputs', variableNames)
    actuator_outputs = actuator_outputs;
elseif ismember('actuator_outputs_1', variableNames)
    actuator_outputs = actuator_outputs_1;
else
    error('Neither actuator_outputs nor actuator_outputs_1 were found in the loaded .mat file.');
end



%% Set the titles
fig = figure;
set(gcf, 'Position', [0, 0, 1250, 600]);

switch selection
    case 1
        exportTitle = 'NED, UVW, $\phi \theta \psi$, PQR';
    case 2
        exportTitle = 'NED, UVW';
    case 3
        exportTitle = 'Phi Theta Psi, PQR';
    case 4
        exportTitle = 'NED';
    case 5
        exportTitle = 'UVW';
    case 6
        exportTitle = 'Phi Theta Psi';
    case 7
        exportTitle = 'PQR';
    case 8
        exportTitle = 'HATCH Quadcopter Controller Performance';
    case 9
        exportTitle = 'HATCH Altitude Hold Performance';
    case 10
        exportTitle = 'RPM Tracking';
    case 11
        exportTitle = '12 DOF Analysis';
    otherwise
        exportTitle = 'Unknown Selection';
end

set(get(gcf,'Children'),'FontSize', 20);


% Set the name of the figure to the name of the .mat file
%set(fig, 'Name', mat_name); % Set the figure name

if selection == 1 % All plots
    num_rows = 4;
    num_cols = 1;
    subplot_pos = [1,1,1,2,2,2,3,3,3,4,4,4];
elseif selection == 2 % just NED and uvw
    num_rows = 2;
    num_cols = 3;
    subplot_pos = [1,2,3,4,5,6,0,0,0,0,0,0];
elseif selection == 3 % just phi,theta,psi, and p,q,r
    num_rows = 2;
    num_cols = 3;
    subplot_pos = [0,0,0,0,0,0,1,2,3,4,5,6];
elseif selection == 4 % just NED
    num_rows = 3;
    num_cols = 1;
    subplot_pos = [1,2,3,0,0,0,0,0,0,0,0,0];
elseif selection == 5 % just u,v,w
    num_rows = 3;
    num_cols = 1;
    subplot_pos = [0,0,0,1,2,3,0,0,0,0,0,0];
elseif selection == 6 % just phi,theta,psi
    num_rows = 3;
    num_cols = 1;
    subplot_pos = [0,0,0,0,0,0,1,2,3,0,0,0];
elseif selection == 7 % just p,q,r
    num_rows = 3;
    num_cols = 1;
    subplot_pos = [0,0,0,0,0,0,0,0,0,1,2,3];
elseif selection == 8 % phi theta r 
    num_rows = 3;
    num_cols = 1;
    subplot_pos = [0,0,0,0,0,0,1,2,0,0,0,3];
elseif selection == 9 % D
    num_rows = 1;
    num_cols = 1;
    subplot_pos = [0,0,1,0,0,0,0,0,0,0,0,0];
elseif selection == 10 % RPM
    num_rows = 1;
    num_cols = 1; 
    subplot_pos = [0,0,0,0,0,0,0,0,0,0,0,0];
elseif selection == 11 % Total plot
    subplot_pos = [0,0,0,0,0,0,0,0,0,0,0,0];
else    
    disp('Invalid Selection')
end

%% x (ft)
if subplot_pos(1) ~= 0
    subplot(num_rows, num_cols, subplot_pos(1));
    h1 = plot(time_seconds, N); hold on;
    ylabel('x (ft)');
    grid on;
end

%% y (ft)

if subplot_pos(2) ~= 0
    subplot(num_rows, num_cols, subplot_pos(2));
    h2 = plot(time_seconds, E);
    ylabel('y (ft)');
    if subplot_pos(4) ~= 0 || subplot_pos(7) ~= 0 || subplot_pos(10) ~= 0
        title('x y z (ft)');
    end
    if selection ~= 4
        xlabel('Time (s)');
    end
    grid on;
end

%% z (ft)
if subplot_pos(3) ~= 0
    subplot(num_rows, num_cols, subplot_pos(3));
    h3 = plot(time_seconds, D);
    if selection == 9 && exist('z_sp_real','var')
        hold on
        
        % Find the indices where ch5_val is 0 (altitude hold mode)
        ch5_val_zero_indices = find(ch5_val == 0);
        
        % Find the indices where ch5_val is 1 (stabilized mode)
        ch5_val_one_indices = find(ch5_val == 1);
        
        % Plot z_sp using dashed lines only for altitude hold mode (ch5_val = 0)
        ch4 = (ch4-1500);
        plot(time_seconds, z_sp_real, 'Color', chicago_maroon);
    
        % Find the indices where the mode changes occur
        mode_change_indices = find(diff(ch5_val) ~= 0);
        
        % Plot vertical lines for mode changes
        for i = 1:length(mode_change_indices)
            if ch5_val(mode_change_indices(i)) == 0
                line_color = 'b';
                mode_label = 'Stabilized';
            elseif ch5_val(mode_change_indices(i)) == 1
                line_color = 'g';
                mode_label = 'Altitude';
            else
                continue; % Skip other modes for now
            end
            
            % Plot vertical line
            plot([time_seconds(mode_change_indices(i)), time_seconds(mode_change_indices(i))], ...
                 [min(min(D,z_sp_)), max(max(D,z_sp_))], line_color, 'LineWidth', 2);
             
            % Label the line
            text(time_seconds(mode_change_indices(i)), min(D), mode_label, ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                 'BackgroundColor', 'w');
        end
        
        legend('$\hat{z}$', '$z_{ref}$')
    end
    if selection == 4
        xlabel('Time (s)')
    end
    ylabel('z (ft)');
    if selection == 1
        ylabel('xyz (ft)')
        xlabel(' ')
        title(' ')
    end
    %ylim([min(min(D)), max(max(D))])
    %xlim([trim_start,end_time])
    grid on;
    
    if selection == 9
        xlabel('Time (s)')
    end

end

if selection == 1 
    for i = 1:length(PTI_switch) 
        hold on
        % Ensure valid range
        if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
            toggleTime = time_seconds(PTI_switch(i));
            xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
            yLimits = ylim(gca);
            % Check the direction of change
            if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                % If rising, label this instance "PTI ON"
                text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            else
                % If falling, label this instance "PTI OFF"
                text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            end
        end
        hold off
    end
    legend([h1,h2,h3],{'$N$','$E$','$D$'})
end

%% u, v, w(ft/s)
% u (ft/s)
if subplot_pos(4) ~= 0
    subplot(num_rows, num_cols, subplot_pos(4));
    h1 = plot(time_seconds, u); hold on;
    ylabel('u (ft/s)');
    grid on;
    
% v (ft/s)
    subplot(num_rows, num_cols, subplot_pos(5));
    h2 = plot(time_seconds, v);
    ylabel('v (ft/s)');
    if subplot_pos(1) ~= 0 || subplot_pos(7) ~= 0 || subplot_pos(10) ~= 0
        title('u v w (ft/s)');
    end
    if selection ~= 5
        xlabel('Time (s)');
    end
    grid on;
    
% w (ft/s)
    subplot(num_rows, num_cols, subplot_pos(6));
    h3 = plot(time_seconds, w);
    xlabel('Time (s)');
    ylabel('w (ft/s)');
    grid on;
    if selection == 1
        ylabel('uvw (ft/s)')
        xlabel(' ')
        title(' ')
    end
end

if selection == 1 
    for i = 1:length(PTI_switch) 
        hold on
        % Ensure valid range
        if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
            toggleTime = time_seconds(PTI_switch(i));
            xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
            yLimits = ylim(gca);
            % Check the direction of change
            if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                % If rising, label this instance "PTI ON"
                text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            else
                % If falling, label this instance "PTI OFF"
                text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            end
        end
        hold off
    end
    legend([h1,h2,h3],{'$u$','$v$','$w$'})
end

%% Phi (deg) 

if subplot_pos(7) ~= 0
    subplot(num_rows, num_cols, subplot_pos(7));
    h1 = plot(time_seconds, phi); hold on;
    if selection == 8 && exist('rollcmd_deg','var')
        plot(time_seconds, rollcmd_deg, 'Color', chicago_maroon);
        for i = 1:length(PTI_switch)
            % Ensure valid range
            if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
                toggleTime = time_seconds(PTI_switch(i));
                xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
                yLimits = ylim(gca);
                % Check the direction of change
                if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                    % If rising, label this instance "PTI ON"
                    text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                else
                    % If falling, label this instance "PTI OFF"
                    text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                end
            end
        end
    end
    if selection ~= 8 && selection ~= 7 && selection ~= 6
        xlabel('Time (s)')
    end
    ylabel('$\phi$ (deg)');
    grid on;
    if selection == 8 && exist('rollcmd_deg','var')
        legend('$\phi$','$\phi_{ref}$','Location','southeast')
        ylim([rollcmd_min-3,rollcmd_max+3])
    end
end

%% Theta (deg)
if subplot_pos(8) ~= 0   
    subplot(num_rows, num_cols, subplot_pos(8));
    h2 = plot(time_seconds, theta); hold on
    if selection == 8 && exist('pitchcmd_deg','var')
        plot(time_seconds, pitchcmd_deg, 'Color', chicago_maroon);
        for i = 1:length(PTI_switch)
            % Ensure valid range
            if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
                toggleTime = time_seconds(PTI_switch(i));
                xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
                yLimits = ylim(gca);
                % Check the direction of change
                if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                    % If rising, label this instance "PTI ON"
                    text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                else
                    % If falling, label this instance "PTI OFF"
                    text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                end
            end
        end
    end
    ylabel('$\theta$ (deg)');
    if selection ~= 8 && selection ~= 7 && selection ~= 6
        xlabel('Time (s)');
    end
    if subplot_pos(1) ~= 0 || subplot_pos(4) ~= 0 || subplot_pos(10) ~= 0
        title('$\phi$ $\theta$ $\psi$ (deg)')
    end
    grid on;
    if selection == 8 && exist('pitchcmd_deg','var')
        legend('$\theta$','$\theta_{ref}$','Location','southeast')
        ylim([pitchcmd_min-3,pitchcmd_max+3])
    end
end

%% Psi (deg)  
if subplot_pos(9) ~= 0
    subplot(num_rows, num_cols, subplot_pos(9));
    h3 = plot(time_seconds, psi);
    xlabel('Time (s)');
    ylabel('$\psi$ (deg)');
    grid on; hold off;
end

if selection == 1 
    for i = 1:length(PTI_switch) 
        hold on
        % Ensure valid range
        if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
            toggleTime = time_seconds(PTI_switch(i));
            xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
            yLimits = ylim(gca);
            % Check the direction of change
            if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                % If rising, label this instance "PTI ON"
                text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            else
                % If falling, label this instance "PTI OFF"
                text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            end
        end
        hold off
    end
    legend([h1,h2,h3],{'$\phi$','$\theta$','$\psi$'})
    if selection == 1
        ylabel('$\phi \theta \psi$ (deg)')
        xlabel(' ')
        title(' ')
    end
end

%% p (deg/s)
if subplot_pos(10) ~= 0
    % Create subplots for p q r
    subplot(num_rows, num_cols, subplot_pos(10));
    h1 = plot(time_seconds, p); hold on;
    if selection ~= 7
        xlabel('Time (s)');
    end
    ylabel('p (deg/s)');
    grid on;
end
    
%% q (deg/s)
if subplot_pos(11) ~= 0
    subplot(num_rows, num_cols, subplot_pos(11));
    h2 = plot(time_seconds, q);
    if selection ~= 7
        xlabel('Time (s)');
    end
    ylabel('q (deg/s)');
    if subplot_pos(1) ~= 0 || subplot_pos(4) ~= 0 || subplot_pos(7) ~= 0
        title('p q r (deg/s)')
    end
    grid on;
end
    
%% r (deg/s)
if subplot_pos(12) ~= 0
    subplot(num_rows, num_cols, subplot_pos(12));
    h3 = plot(time_seconds, r);
    if selection == 8 && exist('yawcmd_dps','var')
        hold on
        plot(time_seconds, yawcmd_dps, 'Color', chicago_maroon);
        for i = 1:length(PTI_switch)
            % Ensure valid range
            if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
                toggleTime = time_seconds(PTI_switch(i));
                xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
                yLimits = ylim(gca);
                % Check the direction of change
                if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                    % If rising, label this instance "PTI ON"
                    text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                else
                    % If falling, label this instance "PTI OFF"
                    text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
                end
            end
        end
        hold off
    end
    xlabel('Time (s)');
    ylabel('r (deg/s)');

    grid on;
    if selection == 8 && exist('yawcmd_dps','var')
        legend('$r$','$r_{ref}$','Location','southeast')
        ylim([yawcmd_min-3,yawcmd_max+3])
    end

end

if selection == 1 
    for i = 1:length(PTI_switch) 
        hold on
        % Ensure valid range
        if PTI_switch(i) > 1 && PTI_switch(i) <= length(ch7_val)
            toggleTime = time_seconds(PTI_switch(i));
            xline(toggleTime, 'b', 'LineWidth', 2); % Plot vertical line
            yLimits = ylim(gca);
            % Check the direction of change
            if ch7_val(PTI_switch(i)+1) > ch7_val(PTI_switch(i))
                % If rising, label this instance "PTI ON"
                text(toggleTime+.1, yLimits(1), 'PTI ON', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            else
                % If falling, label this instance "PTI OFF"
                text(toggleTime+.1, yLimits(1), 'PTI OFF', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Rotation', 90,'FontSize',7,'Color','b');
            end
        end
        hold off
    end
        legend([h1,h2,h3],{'$p$','$q$','$r$'})
    if selection == 1
        ylabel('pqr (deg/s)')
        xlabel(' ')
        title(' ')
    end
end


%% RPM 

if ~exist('esc_status', 'var') && selection == 10
    close all
    error('RPM data is not available in this log file. Please select a different option for the `selection` variable.')
end

if selection == 10 
    % Single motor selection (if applicable)
    motor_number = -1; % -1 - plot all motors | 0 to 3 - plot just that motor

    % Define subplot positions for a 2x2 grid
    subplotPositions = [1, 2, 3, 4];
    
    % Constants for scaling
    OldMin = 376;
    OldMax = 8191;
    NewMin = 500.3540471249;
    NewMax = 10900;
    RPM_drop_threshold = NewMin;
    
    % Loop through each motor, now accounting for 0:3 indexing
    for motorIndex = 0:3
        if motor_number >= 0 && motorIndex ~= motor_number
            continue;
        end
        if motor_number == -1
            subplot(2, 2, subplotPositions(motorIndex + 1)); % Adjusted for 0-based indexing
        else
            subplot(1,1,subplotPositions(motorIndex + 1)); % Adjusted for 0-based indexing
        end
        
        % Initialize variables to store RPM commands from different sources
        RPM_cmd_motor_UAVTlbx = []; % For UAV Toolbox method
        RPM_cmd_motor_ActuatorOutputs = []; % For ActuatorOutputs method

        % Check if 'rpm_cmd' is available in 'custom_internal'
        if exist('custom_internal','var') && istimetable(custom_internal) && ismember('rpm_cmd', custom_internal.Properties.VariableNames)
            % Accessing commanded RPM for the current motor via custom_internal
            RPM_cmd_motor_UAVTlbx = custom_internal.rpm_cmd(:, motorIndex+1); % Adjusted for 0-based indexing
        end
        
        if ~exist('custom_internal','var') && istimetable(custom_internal) && ismember('rpm_cmd', custom_internal.Properties.VariableNames)
            % Accessing commanded RPM for the current motor via fallback method
            RPM_cmd_motor_ActuatorOutputs = eval(sprintf('actuator_outputs.output(:,%d)', motorIndex + 1)); % Adjusted for 0-based indexing, MATLAB indexing starts at 1
            RPM_cmd_motor_ActuatorOutputs(RPM_cmd_motor_ActuatorOutputs == 65535) = 0; % Handle specific cases
    
            % Scale the ActuatorOutputs RPM values
            RPM_cmd_motor_ActuatorOutputs = (RPM_cmd_motor_ActuatorOutputs - OldMin) * (NewMax - NewMin) / (OldMax - OldMin) + NewMin;
            
            % Trim and handle data for both commanded methods based on start and end index
            RPM_cmd_motor_ActuatorOutputs_trimmed = RPM_cmd_motor_ActuatorOutputs(start_index:end_index);
        end
        
        
        % Accessing actual RPM for the current motor
        RPM_actual_motor = eval(sprintf('esc_status.("esc[%d].esc_rpm")', motorIndex)); % Adjusted for 0-based indexing
        % Trim actual motor data
        RPM_actual_motor_trimmed = RPM_actual_motor(start_index:end_index);

        if ~isempty(RPM_cmd_motor_UAVTlbx)
            RPM_cmd_motor_UAVTlbx_trimmed = RPM_cmd_motor_UAVTlbx(start_index:end_index);
            RPM_cmd_motor_UAVTlbx_trimmed(RPM_cmd_motor_UAVTlbx_trimmed == 65535) = 0; % Handle specific cases
        end
                
        % Plotting
        if ~isempty(RPM_cmd_motor_UAVTlbx)
            plot(time_seconds, RPM_cmd_motor_UAVTlbx_trimmed, 'Color', burnt_orange, 'LineWidth', 1, 'DisplayName', 'RPM cmd + SysID Signal'); % Green, thicker dotted line for UAV Toolbox Command
            hold on;
        end

        if isempty(RPM_cmd_motor_UAVTlbx)
            plot(time_seconds, RPM_cmd_motor_ActuatorOutputs_trimmed, '--', 'Color', burnt_orange, 'LineWidth', 1, 'DisplayName', 'ActuatorOutputs Cmd'); % Thicker dashed line
            hold on;
        end
        sysIDInput = custom_internal.sysid_uout(:,motorIndex+1);
        sysIDTotal = sysIDInput(start_index:end_index)-RPM_cmd_motor_UAVTlbx_trimmed;
        %plot(time_seconds, sysIDTotal,'--', 'Color', 'b', 'LineWidth',1,'DisplayName', 'RPM cmd')
        plot(time_seconds, RPM_actual_motor_trimmed, '-.','Color', chicago_maroon, 'LineWidth', 1, 'DisplayName', 'RPM Actual');
        ylim([(0) (12500)])
        xlim([trim_start,end_time])
        title(sprintf('Motor %d', motorIndex+1));
        xlabel('Time (s)');
        ylabel('RPM');
        grid on;
        legend('Location', 'northwest','FontSize',10)
    end
end

%% Total figure
time_seconds = seconds(time_seconds); % Transpose to column vector

if selection == 11
    % Ensure the required data is available
    if exist('vehicle_local_position', 'var') && ...
       exist('vehicle_angular_velocity', 'var') && ...
       exist('vehicle_attitude', 'var')

        % Create timetables directly using column vectors
        stateData = timetable(time_seconds, [N,E,D],[phi,theta,psi],[u,v,w],[p,q,r], 'VariableNames', {'NED_ft','EulerAngles_deg','uvw_ft_s','omega_deg_s'});

        % Create figure for stacked plots
        fig;
        stackedplot(stateData,{'NED_ft','EulerAngles_deg','uvw_ft_s','omega_deg_s'})
        
    else
        error('One or more required variables (vehicle_local_position, vehicle_angular_velocity, vehicle_attitude) are missing.');
    end
end

%% Exporting
if export == true
    exportgraphics(fig,[exportTitle,'.pdf'],'ContentType','vector')
end

end