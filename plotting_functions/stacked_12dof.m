function fig = stacked_12dof(time_seconds, all_12dof)
% 12 DOF Analysis Visualization
%
% DESCRIPTION:
%   This function visualizes a 12 degree of freedom analysis, including 
%   position, orientation, velocity, and angular velocity over time.
%
% INPUTS:
%   time_seconds - Time vector in seconds
%   all_12dof    - 12 DOF data matrix [num_samples x 12]
%                  Columns: [N, E, D, phi, theta, psi, u, v, w, p, q, r]
%
% OUTPUTS:
%   Figure displaying the 12 DOF analysis
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

% Extract 12 DOF data
N = all_12dof(:,1); 
E = all_12dof(:,2); 
D = all_12dof(:,3);
phi = all_12dof(:,4); 
theta = all_12dof(:,5); 
psi = all_12dof(:,6);
u = all_12dof(:,7); 
v = all_12dof(:,8); 
w = all_12dof(:,9);
p = all_12dof(:,10); 
q = all_12dof(:,11); 
r = all_12dof(:,12);

% Create timetables directly using column vectors
stateData = timetable(time_seconds, [N, E, D], [phi, theta, psi], [u, v, w], [p, q, r], ...
                      'VariableNames', {'NED_ft', 'EulerAngles_deg', 'uvw_ft_s', 'omega_deg_s'});

% Create figure with stacked plots
stackedplot(stateData, {'NED_ft', 'EulerAngles_deg', 'uvw_ft_s', 'omega_deg_s'});

end