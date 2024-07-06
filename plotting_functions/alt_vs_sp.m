function fig = alt_vs_sp(time_seconds, alt_ft, alt_ft_sp)

% Altitude vs Altitude Setpoint Visualization
%
% DESCRIPTION:
%   This function visualizes the altitude and the altitude setpoint over time.
%
% INPUTS:
%   time_seconds - Time vector in seconds
%   alt_ft       - Altitude in feet
%   alt_ft_sp    - Altitude setpoint in feet
%
% OUTPUTS:
%   Figure displaying the altitude and altitude setpoint over time
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

plot(time_seconds, alt_ft, 'DisplayName', 'Altitude');
hold on;
plot(time_seconds, alt_ft_sp, 'DisplayName', 'Altitude Setpoint');

xlabel('Time (s)');
ylabel('Altitude (ft)');
legend('Location', 'northwest');
grid on;

end