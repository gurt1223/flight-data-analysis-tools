function log = trimLogData(log, start_idx, end_idx, time_seconds)

% Trim Log Data to Specified Indices
%
% DESCRIPTION:
%   This function dynamically trims all fields in the log structure based
%   on specified start and end indices. It handles special cases for
%   quaternion and vertical velocity setpoint fields.
%
% INPUTS:
%   log          - Structure containing various log data fields
%   start_idx    - Start index for trimming
%   end_idx      - End index for trimming
%   time_seconds - Time vector in seconds
%
% OUTPUTS:
%   log          - Trimmed log structure with updated fields
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

    % List of fields that need special handling
    specialFields = {'quat', 'vv_fps_sp'};
    
    % Loop through all fields in the log structure
    fields = fieldnames(log);
    for i = 1:numel(fields)
        field = fields{i};
        
        if ismember(field, specialFields)
            % Special handling for specific fields
            switch field
                case 'quat'
                    log.(field) = log.(field)(start_idx:end_idx, :);
                case 'vv_fps_sp'
                    log.alt_ft_sp = cumtrapz(time_seconds, log.vv_fps_sp(start_idx:end_idx));
                    log.vv_fps_sp = log.vv_fps_sp(start_idx:end_idx);
            end
        else
            % General handling for all other fields
            if size(log.(field), 1) >= end_idx
                if size(log.(field), 1) == size(time_seconds, 1)
                    log.(field) = log.(field)(start_idx:end_idx, :);
                else
                    log.(field) = log.(field)(start_idx:end_idx,:);
                end
            end
        end
    end
    
    % Convert to Euler Angles (deg)
    eulXYZ = quat2eul(log.quat, "ZYX");
    log.phi = eulXYZ(:, 3) * (180 / pi);
    log.theta = eulXYZ(:, 2) * (180 / pi);
    log.psi = eulXYZ(:, 1) * (180 / pi);

    % Make large 12dof variable
    log.all_12dof = [log.N, log.E, log.D, log.phi, log.theta, log.psi, log.u, log.v, log.w, log.p, log.q, log.r];
end
