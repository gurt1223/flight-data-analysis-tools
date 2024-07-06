function log = createLogStructure(mat_name, mat_path)

% Create Log Structure from MAT File
%
% DESCRIPTION:
%   This function loads data from a specified .mat file and assembles it
%   into a structured format. The log structure includes various fields
%   such as position, velocity, quaternion, angular velocity, vertical
%   velocity, altitude, and innerloop commands. The user must ensure that
%   the appropriate parts of the ulog file are linked to the specified
%   fields so that the main plotting script functions correctly.
%
% INPUTS:
%   mat_name  - Name of the .mat file containing the flight data
%   mat_path  - Path to the directory containing the .mat file
%
% OUTPUTS:
%   log       - Structure containing the flight data with specified fields
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

    % Load the .mat file
    load(fullfile(mat_path, mat_name));
    
    m2ft = 3.28084; % meters to feet conversion

    % Position (ft)
    log.N = (vehicle_local_position.x)*m2ft;
    log.E = (vehicle_local_position.y)*m2ft;
    log.D = (vehicle_local_position.z)*m2ft;

    % Velocity (ft/s)
    log.u = (vehicle_local_position.vx)*m2ft;
    log.v = (vehicle_local_position.vy)*m2ft;
    log.w = (vehicle_local_position.vz)*m2ft;

    % Quaternion (quat = [q0,q1,q2,q3])
    log.quat = vehicle_attitude.q;

    % Angular Velocity (deg/s)
    log.p = vehicle_angular_velocity.xyz(:,1)*(180/pi);
    log.q = vehicle_angular_velocity.xyz(:,2)*(180/pi);
    log.r = vehicle_angular_velocity.xyz(:,3)*(180/pi);

    % Vertical Velocity (ft/s)
    log.vv_fps = fcs_signals.px4_data(:,9); 

    % Altitude (ft)
    log.alt_ft = fcs_signals.px4_data(:,13);

    % Vertical Velocity Setpoint (ft/s)
    log.vv_fps_sp = fcs_signals_aux.fcs_data(:,16); 

    % Innerloop commands (deg || deg/s)
    log.phi_sp = fcs_signals.att_cmd(:,1)*(180/pi); % roll angle (deg)
    log.theta_sp = fcs_signals.att_cmd(:,2)*(180/pi); % pitch angle (deg)
    log.r_sp = fcs_signals.att_cmd(:,3)*(180/pi); % yaw angular rate (deg/s)

    num_motors = 4; % Number of motors / signals present on esc_status

    % Preallocate the log.rpm and log.rpm_sp arrays
    log.rpm = zeros(size(esc_status.('esc[0].esc_rpm'), 1), num_motors);
    log.rpm_sp = zeros(size(fcs_signals.eng_cmd, 1), num_motors);
    
    for motorIndex = 1:num_motors
        log.rpm(:, motorIndex) = esc_status.(sprintf('esc[%d].esc_rpm', motorIndex));
        log.rpm_sp(:, motorIndex) = fcs_signals.eng_cmd(:, motorIndex);
    end
end
