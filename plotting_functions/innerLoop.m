function fig = innerLoop(time_seconds, inner, inner_sp)
    
% Inner Loop Performance Visualization
%
% DESCRIPTION:
%   This function visualizes the roll angle, pitch angle, and yaw rate 
%   along with their commanded values over time.
%
% INPUTS:
%   time_seconds - Time vector in seconds
%   inner        - Actual roll angle, pitch angle, and yaw rate [phi, theta, r]
%   inner_sp     - Commanded roll angle, pitch angle, and yaw rate [phi_sp, theta_sp, r_sp]
%
% OUTPUTS:
%   Figure displaying the inner loop performance
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

% Extract the inner loop data
phi = inner(:, 1);
theta = inner(:, 2);
r = inner(:, 3);

% Extract the setpoint data
phi_sp = inner_sp(:, 1);
theta_sp = inner_sp(:, 2);
r_sp = inner_sp(:, 3);

% Plot roll angle (phi)
subplot(3,1,1)
plot(time_seconds, phi, 'DisplayName', '$\phi$'); hold on;
plot(time_seconds, phi_sp, 'DisplayName', '$\phi_{ref}$');
ylabel('$\phi$ (deg)', 'Interpreter', 'latex');
grid on;
legend('Location', 'southeast', 'Interpreter', 'latex');

% Plot pitch angle (theta)
subplot(3,1,2)
plot(time_seconds, theta, 'DisplayName', '$\theta$'); hold on;
plot(time_seconds, theta_sp, 'DisplayName', '$\theta_{ref}$');
ylabel('$\theta$ (deg)', 'Interpreter', 'latex');
grid on;
legend('Location', 'southeast', 'Interpreter', 'latex');

% Plot yaw rate (r)
subplot(3,1,3)
plot(time_seconds, r, 'DisplayName', '$r$'); hold on;
plot(time_seconds, r_sp, 'DisplayName', '$r_{ref}$');
xlabel('Time (s)');
ylabel('r (deg/s)');
grid on;
legend('Location', 'southeast', 'Interpreter', 'latex');

end
