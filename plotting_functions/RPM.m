function fig = RPM(time_seconds, rpm, rpm_sp)

% RPM Tracking Visualization
%
% DESCRIPTION:
%   This function visualizes the RPM and the RPM setpoint for multiple motors over time.
%
% INPUTS:
%   time_seconds - Time vector in seconds
%   num_motors   - Number of motors
%   rpm          - Actual RPM values [num_samples x num_motors]
%   rpm_sp       - Commanded RPM values [num_samples x num_motors]
%
% OUTPUTS:
%   Figure displaying the RPM tracking for multiple motors
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

% Define the figure 
fig = figure;

% Single motor selection (if applicable)
motor_number = -1; % -1: plot all motors | 0 to 3: plot specific motor

% Calculate the number of motors
num_motors = size(rpm, 2);

% Calculate the number of rows and columns for the subplot grid
num_cols = ceil(sqrt(num_motors));
num_rows = ceil(num_motors / num_cols);

% Loop through each motor
for motor_Idx = 1:num_motors
    if motor_number >= 0 && motor_Idx ~= motor_number
        continue;
    end
    
    % Create a subplot for each motor
    if motor_number == -1
        subplot(num_rows, num_cols, motor_Idx); 
    else
        subplot(1, 1, motor_Idx); 
    end

    plot(time_seconds, rpm_sp(:,motor_Idx), '-.', 'DisplayName', 'RPM Setpoint');
    hold on;
    plot(time_seconds, rpm(:,motor_Idx), 'DisplayName', 'RPM Actual');
    hold off;
    
    xlim([time_seconds(1), time_seconds(end)]);
    title(sprintf('Motor %d', motor_Idx));
    xlabel('Time (s)');
    ylabel('RPM');
    grid on;
    legend('Location', 'southwest', 'FontSize', 10);
end

end
