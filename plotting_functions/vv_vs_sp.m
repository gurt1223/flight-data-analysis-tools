function fig = vv_vs_sp(time_seconds, vv_fps, vv_fps_sp)

% Vertical Velocity vs Setpoint Visualization
%
% DESCRIPTION:
%   This function visualizes the vertical velocity and its setpoint over time.
%
% INPUTS:
%   time_seconds - Time vector in seconds
%   vv_fps       - Vertical velocity in feet per second
%   vv_fps_sp    - Vertical velocity setpoint in feet per second
%
% OUTPUTS:
%   Figure displaying the vertical velocity vs setpoint
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

plot(time_seconds, vv_fps, 'DisplayName', 'Vertical Velocity');
hold on;
plot(time_seconds, vv_fps_sp, 'DisplayName', 'Vertical Velocity Setpoint');
hold off;

xlabel('Time (s)');
ylabel('Vertical Velocity (ft/s)');
legend('Location', 'best');
grid on;

end