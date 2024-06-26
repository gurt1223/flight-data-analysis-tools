function processUlogTimeSynced(ulog_path, ulog_name, freq, CutoffFrequency)
% Convert ULOG to MAT file and time synchronize the data
%
% (C) 2024 Garrett D. Asper <garrettasper@vt.edu>
% A Ulog file for analysis
%
% Credit: Filtering elements of this code are thanks to Jeremy Hopwood
% from jwhgit01 FlightTest-Tools/src/FlightDataProcessing/FlightData.m
% 
% Description: This script can be run in two modes:
%   - Interactive Mode: Without input arguments, it prompts the user to select
%     a ULOG file for conversion through a graphical interface.
%   - Programmatic Mode: With provided ulog_path and ulog_name arguments, it directly
%     processes the specified ULOG file without user interaction.
%
% MATLAB Requirements: Matlab R2023a or newer, UAV Toolbox
%
% The following topics will be time-synchronized if available in the log
% file:
%
%       actuator_controls_0
%       actuator_outputs
%       actuator_motors
%       battery_status
%       cpuload
%       custom_internal
%       esc_status
%       input_rc
%       sensor_combined
%       sensor_gps
%       sensor_rpm
%       vehicle_acceleration
%       vehicle_air_data
%       vehicle_angular_velocity
%       vehicle_attitude
%       vehicle_attitude_setpoint
%       vehicle_control_mode
%       vehicle_global_position
%       vehicle_local_position
%       vehicle_local_position_setpoint
%       vehicle_rates_setpoint
%       vehicle_status
%       vehicle_torque_setpoint
%       vehicle_thrust_setpoint
%
% Input(s): 
%     ulog_path (optional): Path to the directory containing the .ulg file
%     ulog_name (optional): Name of the .ulg file to be processed
%     freq (optional): Frequency for interpolation
%
% Output(s): .mat file containing time-synchronized data
%
% Usage:
%   Interactive Mode: processUlogTimeSynced
%   Programmatic Mode: processUlogTimeSynced('C:/data', 'flightData.ulg')

% Error checks
%% Input Validation and Error Handling

% ~~~~ Interactive Mode ~~~~ (Define these and then run the script)
    if nargin == 0 || isempty(ulog_name) || isempty(ulog_path) 
        clear; clc; % Use clear if you want to reset the environment in interactive mode
        % Prompt the user to select the Ulog file using the file explorer
        ulog_path = strrep(mfilename('fullpath'), mfilename, '');
        [ulog_name, ulog_path] = uigetfile([ulog_path, '*.ulg'], 'MultiSelect', 'off');
        ulog_name = {ulog_name};
        CutoffFrequency = 6; % [Hz]
        ESC_filter = false; % determines whether ESC signals will be filtered

        % Initialize parameters
        freq = 200; % Desired data frequency in Hz
    end

% ~~~~ Programatic Mode ~~~~ (This will run if function is called
% externally)
if nargin == 4
    ulog_name = {ulog_name};
    ESC_filter = false;
end

tic % starts the timer

% Validate the provided or selected file
if isequal(ulog_name,0) || isequal(ulog_path,0)
    error('File selection was cancelled or invalid file was provided.');
end

dt = 1 / freq; % Time step

% Import Ulog
ulog_obj = ulogreader([ulog_path, ulog_name{1}]);

% Read Ulog data
ulogData = readTopicMsgs(ulog_obj);

% Extract topic data
dataTopicNames = ulogData.TopicNames;
disp(' ')
disp('************************************************************')
disp('The following topics are present on the log file:')
disp(' ')
disp(dataTopicNames)
dataTopicMessages = ulogData.TopicMessages;
PX4_parameters = readParameters(ulog_obj);

% Prepare to list uninterpolated topics
disp(' ')
disp('************************************************************')
disp('The following topics were present but have not been retimed:')
disp(' ')
    
% Process and store data as individual variables
variableNames = cell(size(dataTopicNames));

% Initialize vector to store topic sample rates
fs_vec = [0];
topicNames_vec = {};

for i = 1:length(dataTopicNames)
    % Handle multiple topics with the same name
    topicName = dataTopicNames{i};
    indices = find(strcmp(topicName, dataTopicNames));
    if length(indices) == 1
        extra = '';
    elseif length(indices) > 1
        extra = ['_' num2str(find(indices == i))];
    else
        error('Data field does not exist');
    end

    variableNames{i} = [topicName, extra];
    eval([variableNames{i}, ' = dataTopicMessages{i};']) % Store data as timetable
end

% Define the structure that all retimed data will be saved to
retimedData = struct();

    % Retiming section of the code
    % logger start time
    t0_ULOG = seconds(ulog_obj.StartTime);
    tf_ULOG = seconds(ulog_obj.EndTime);
    
    % Define time
    Time = seconds(t0_ULOG:dt:tf_ULOG).';
   
    % Perform retiming if the topic is present

    % actuator_controls_0
    if any(strcmp(ulogData.TopicNames,'actuator_controls_0'))
        fs_vec(end+1) = 1/median(diff(seconds(actuator_controls_0.timestamp)));
        topicNames_vec{end+1} = 'actuator_controls_0';
        actuator_controls_0_retimed = retime(actuator_controls_0, Time, 'pchip');
        actuator_controls_0 = actuator_controls_0_retimed;
    end

    % actuator_motors
    if any(strcmp(ulogData.TopicNames,'actuator_motors'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(actuator_motors.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'actuator_motors';

        idx = find(strcmp(ulogData.TopicNames, 'actuator_motors'));

        % Access the corresponding TopicMessages using the index
        actuator_motors_data = ulogData.TopicMessages{idx};
        % Extract the relevant data (assuming 'states' and 'covariances' are table variables)
        
        actuator_motors = actuator_motors_data(:,{'control'});

        actuator_motors_retimed = retime(actuator_motors, Time, 'pchip');
        clear actuator_motors_data
        actuator_motors = actuator_motors_retimed;
    end

    % actuator_outputs
    if any(strcmp(ulogData.TopicNames,'actuator_outputs'))
        if exist('actuator_outputs', 'var') == 1
            % Find the original sample rate
            fs_vec(end+1) = 1/median(diff(seconds(actuator_outputs.timestamp))); % Sampling frequency
            topicNames_vec{end+1} = 'actuator_outputs';

            % Retime
            actuator_outputs_retimed = retime(actuator_outputs, Time, 'nearest');
            actuator_outputs = actuator_outputs_retimed;
        end

        if exist('actuator_outputs_1', 'var') == 1
            % Find the original sample rate
            fs_vec(end+1) = 1/median(diff(seconds(actuator_outputs_1.timestamp))); % Sampling frequency
            topicNames_vec{end+1} = 'actuator_outputs_1';

            % Retime
            actuator_outputs_retimed = retime(actuator_outputs_1, Time, 'nearest');
            actuator_outputs_1 = actuator_outputs_retimed;
        end

        if exist('actuator_outputs_2', 'var') == 1
            % Find the original sample rate
            fs_vec(end+1) = 1/median(diff(seconds(actuator_outputs_2.timestamp))); % Sampling frequency
            topicNames_vec{end+1} = 'actuator_outputs_2';
            
            % Retime
            actuator_outputs_retimed = retime(actuator_outputs_2, Time, 'nearest');
            actuator_outputs_2 = actuator_outputs_retimed;
        end
        
    end
    
    % battery_status
    if any(strcmp(ulogData.TopicNames,'battery_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(battery_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'battery_status';

        % Retime
        battery_status_retimed = retime(battery_status, Time, 'nearest');
        battery_status = battery_status_retimed;
    end
    
    % cpuload
    if any(strcmp(ulogData.TopicNames,'cpuload'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(cpuload.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'cpuload';

        % Retime
        cpuload_retimed = retime(cpuload, Time, 'pchip');
        cpuload = cpuload_retimed;
    end
    
    % custom_internal
    if any(strcmp(ulogData.TopicNames,'custom_internal'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(custom_internal.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'custom_internal';

        % Retime
        custom_internal = convertvars(custom_internal,vartype('numeric'),'double');
        custom_internal_retimed = retime(custom_internal, Time, 'pchip');
        custom_internal = custom_internal_retimed;
    end

    % esc_status
    if any(strcmp(ulogData.TopicNames,'esc_status'))
        % Find the index of 'esc_status' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'esc_status'));
    
        % Access the corresponding TopicMessages using the index
        esc_status = ulogData.TopicMessages{idx};
        timestamp_ = seconds(seconds(esc_status.timestamp));
        
        % Convert all variables in esc_status to double
        allVars = vartype('numeric'); % Selects all numeric variables
        esc_status = convertvars(esc_status, allVars, 'double');
        
        % Get number of ESCs
        esc_count = esc_status.esc_count(1,1);
    
        for ii = 1:esc_count
            % Process RPM and timestamp for each ESC as your colleague did
            ni = double(esc_status.(['esc[' num2str(ii-1) '].esc_rpm']));
            % Remove duplicate timestamp values and keep first sample
            [ti,ia,~] = unique(esc_status.(['esc[' num2str(ii-1) '].timestamp']),'first');
            tidd = double(ti)/1e6;
            nidd = ni(ia);
            esc_rpm = nidd; % defined the current ESCs RPMs over time
            % Calculate the sampling frequency
            timeInSeconds = tidd; % Ensure Time is in seconds
            fs = 1/median(diff(timeInSeconds)); % Sampling frequency
            fs_esc_status(ii) = fs;
            
            % Check if filtering is needed based on cutoff frequency and sampling rate
            if floor(fs) > 2*CutoffFrequency && ESC_filter == true
                Wn = 2*CutoffFrequency/fs; % Fraction of the Nyquist rate
                [b, a] = butter(5, Wn); % 5th order low-pass Butterworth filter

                % Apply the filter to each column of esc_rpm
                esc_rpm = filtfilt(b, a, nidd);
            end
        
            % Finally, assign retimed and possibly filtered RPM data back to esc_status
            esc_rpm_ = timetable(seconds(tidd), esc_rpm);
            esc_status_retimed(:,ii) = retime(esc_rpm_,Time, 'pchip');
        end
    
            % Generate new names for the esc_rpm columns
            numMotors = 4; % Adjust this based on how many 'esc_rpmX' columns you have
            for motorIndex = 0:(numMotors-1)
                newVariableNames{motorIndex+1} = sprintf('esc[%d].esc_rpm', motorIndex);
            end
           
            % Apply the new variable names to the timetable
            esc_status = renamevars(esc_status_retimed, allVars, newVariableNames);
            fs_vec(end+1) = median(fs_esc_status);
            topicNames_vec{end+1} = 'esc_status';
    end


    % estimator_states
    if any(strcmp(ulogData.TopicNames, 'estimator_states'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(estimator_states.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'estimator_states';
   
        % Find the index of 'estimator_states' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'estimator_states'));
    
        % Access the corresponding TopicMessages using the index
        estimator_states_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data (assuming 'states' and 'covariances' are table variables)
        estimator_states = estimator_states_data(:, {'states', 'covariances'});
    
        % Perform retiming
        estimator_states_retimed = retime(estimator_states, Time, 'pchip');
        clear estimator_states_data
        estimator_states = estimator_states_retimed;
    end
    
    % estimator_status
    if any(strcmp(ulogData.TopicNames,'estimator_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(estimator_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'estimator_status';

        % Find the index of 'estimator_status' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'estimator_status'));
    
        % Access the corresponding TopicMessages using the index
        estimator_status_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data (assuming 'states' and 'covariances' are table variables)
        estimator_status = estimator_status_data(:,...
                {'pos_horiz_accuracy','pos_vert_accuracy'});
        estimator_status_retimed = retime(estimator_status, Time, 'pchip');
        clear estimator_status_data
        estimator_status = estimator_status_retimed;
    end
    
    % estimator_event_flags
    if any(strcmp(ulogData.TopicNames,'estimator_event_flags_1')) || any(strcmp(ulogData.TopicNames,'estimator_event_flags'))
        if any(strcmp(ulogData.TopicNames,'estimator_event_flags_1'))
            % Find the original sample rate
            fs_vec(end+1) = 1/median(diff(seconds(estimator_event_flags_1.timestamp))); % Sampling frequency
            topicNames_vec{end+1} = 'estimator_event_flags_1';

            % Retime
            estimator_event_flags_retimed = retime(estimator_event_flags_1, Time, 'nearest');
        else
            % Find the original sample rate
            fs_vec(end+1) = 1/median(diff(seconds(estimator_event_flags.timestamp))); % Sampling frequency
            topicNames_vec{end+1} = 'estimator_event_flags';

            % Retime
            estimator_event_flags_retimed = retime(estimator_event_flags, Time, 'nearest');
        end
        estimator_event_flags = estimator_event_flags_retimed;
    end

    % event
    if any(strcmp(ulogData.TopicNames,'event'))        
        % Retiming has not yet been incorporated for this topic
        disp('event')
    end
    
    % input_rc
    if any(strcmp(ulogData.TopicNames,'input_rc'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(input_rc.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'input_rc';
    
        % Retime
        input_rc_retimed = retime(input_rc, Time, 'nearest');
        input_rc = input_rc_retimed;
    end
        
    % sensor_combined
    if any(strcmp(ulogData.TopicNames,'sensor_combined'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(sensor_combined.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'sensor_combined';

        % Find the index of 'sensor_combined' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'sensor_combined'));
    
        % Access the corresponding TopicMessages using the index
        sensor_combined_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data (assuming 'gyro_rad' and 'accelerometer_m_s2' are variables)
        sensor_combined = sensor_combined_data(:,...
            {'gyro_rad','accelerometer_m_s2'});
    
        % Calculate sampling frequency
        fs = 1 / median(diff(seconds(sensor_combined.timestamp)));
    
        % Filter data if applicable
        if ~isempty(CutoffFrequency) && (floor(fs) > 2 * CutoffFrequency)
            Wn = 2 * CutoffFrequency / fs; % fraction of the Nyquist rate
            [b, a] = butter(5, Wn); % 5th order lowpass Butterworth filter
            sensor_combined.gyro_rad = filtfilt(b, a, sensor_combined.gyro_rad);
            sensor_combined.accelerometer_m_s2 = filtfilt(b, a, sensor_combined.accelerometer_m_s2);
        end
    
        % Re-time the data
        sensor_combined_retimed = retime(sensor_combined, Time, 'pchip', 'EndValues', 'extrap');
        sensor_combined = sensor_combined_retimed;
        clear sensor_combined_data
    end

    
    % sensor_gps
    if any(strcmp(ulogData.TopicNames,'sensor_gps'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(sensor_gps.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'sensor_gps';

        % Find the index of 'sensor_gps' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'sensor_gps'));
    
        % Access the corresponding TopicMessages using the index
        sensor_gps_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data (assuming 'states' and 'covariances' are table variables)
        sensor_gps = sensor_gps_data(:,{'time_utc_usec',...
        'lat','lon','alt','alt_ellipsoid','s_variance_m_s','eph','epv',...
        'vel_n_m_s','vel_e_m_s','vel_d_m_s'});
        sensor_gps = removevars(sensor_gps,'time_utc_usec');
        % Convert integers to appropriate units
        sensor_gps.lat = double(sensor_gps.lat)*1e-7;
        sensor_gps.lon = double(sensor_gps.lon)*1e-7;
        sensor_gps.alt = double(sensor_gps.alt)*1e-3;
        sensor_gps.alt_ellipsoid = double(sensor_gps.alt_ellipsoid)*1e-3;
        sensor_gps_retimed = retime(sensor_gps, Time, 'pchip');
        clear sensor_gps_data
        sensor_gps = sensor_gps_retimed;
    end
    
    % sensor_rpm
    if any(strcmp(ulogData.TopicNames,'sensor_rpm'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(sensor_rpm.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'sensor_rpm';
        
        % Retime
        sensor_rpm_retimed = retime(sensor_rpm, Time, 'nearest');
        disp('sensor_rpm needs revisited if needed')
    end

    % vehicle_acceleration
    if any(strcmp(ulogData.TopicNames,'vehicle_acceleration'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_acceleration.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_acceleration';

        % Retime
        vehicle_acceleration_retimed = retime(vehicle_acceleration, Time, 'pchip');
        vehicle_acceleration = vehicle_acceleration_retimed;
    end
    
    % vehicle_air_data
    if any(strcmp(ulogData.TopicNames,'vehicle_air_data'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_air_data.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_air_data';

        disp('vehicle_air_data')
    end
    
    % vehicle_angular_velocity
    if any(strcmp(ulogData.TopicNames,'vehicle_angular_velocity'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_angular_velocity.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_angular_velocity';

        % Retime
        vehicle_angular_velocity_retimed = retime(vehicle_angular_velocity, Time, 'pchip');
        vehicle_angular_velocity = vehicle_angular_velocity_retimed;
    end
    
    % vehicle_attitude
    if any(strcmp(ulogData.TopicNames,'vehicle_attitude'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_attitude.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_attitude';

        % Find the index of 'vehicle_attitude' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'vehicle_attitude'));
    
        % Access the corresponding TopicMessages using the index
        vehicle_attitude_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data 
        vehicle_attitude = vehicle_attitude_data(:,'q');
        vehicle_attitude_retimed = retime(vehicle_attitude, Time, 'pchip');
        clear vehicle_attitude_Data
        vehicle_attitude = vehicle_attitude_retimed;
    end
    
    % vehicle_attitude_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_attitude_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_attitude_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_attitude_setpoint';

        % Retime
        vehicle_attitude_setpoint_retimed = retime(vehicle_attitude_setpoint, Time, 'previous');
        vehicle_attitude_setpoint = vehicle_attitude_setpoint_retimed;
    end
    
    % vehicle_control_mode
    if any(strcmp(ulogData.TopicNames,'vehicle_control_mode'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_control_mode.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_control_mode';

        % Retime
        vehicle_control_mode_retimed = retime(vehicle_control_mode, Time, 'previous');
        vehicle_control_mode = vehicle_control_mode_retimed;
    end
    
    % vehicle_global_position
    if any(strcmp(ulogData.TopicNames,'vehicle_global_position'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_global_position.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_global_position';

        % Retime
        vehicle_global_position = removevars(vehicle_global_position,'lat_lon_reset_counter');
        vehicle_global_position = removevars(vehicle_global_position,'alt_reset_counter');
        vehicle_global_position = removevars(vehicle_global_position,'terrain_alt_valid');
        vehicle_global_position = removevars(vehicle_global_position,'dead_reckoning');
        vehicle_global_position_retimed = retime(vehicle_global_position, Time, 'pchip');
        vehicle_global_position = vehicle_global_position_retimed;
    end
    
    % vehicle_local_position
    if any(strcmp(ulogData.TopicNames,'vehicle_local_position'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_local_position.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_local_position';

        % Find the index of 'vehicle_local_position' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'vehicle_local_position'));
    
        % Access the corresponding TopicMessages using the index
        vehicle_local_position_data = ulogData.TopicMessages{idx};
    
        % Extract the relevant data 
        vehicle_local_position = vehicle_local_position_data(:,{'x','y','z',...
        'vx','vy','vz','ax','ay','az','heading'});
        vehicle_local_position_retimed = retime(vehicle_local_position, Time, 'pchip');
        clear vehicle_local_position_Data
        vehicle_local_position = vehicle_local_position_retimed;
    end
    
    % vehicle_local_position_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_local_position_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_local_position_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_local_position_setpoint';
        
        % Retime
        vehicle_local_position_setpoint_retimed = retime(vehicle_local_position_setpoint(:,["x","y","z","vx","vy","vz"]),Time,'previous');
        vehicle_local_position_setpoint = vehicle_local_position_setpoint_retimed;
    end
    
    % vehicle_rates_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_rates_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_attitude.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_rates_setpoint';

        disp('vehicle_rates_setpoint')
    end
    
    % vehicle_status
    if any(strcmp(ulogData.TopicNames,'vehicle_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_status';
        
        % Retime
        vehicle_status_retimed = retime(vehicle_status, Time, 'nearest');
        vehicle_status = vehicle_status_retimed;
    end
    
    % vehicle_torque_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_torque_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_torque_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_torque_setpoint';
        disp('vehicle_torque_setpoint')
    end
    
    % vehicle_thrust_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_thrust_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_thrust_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_thrust_setpoint';
        disp('vehicle_thrust_setpoint')
    end
    
    % yaw_estimator_status
    if any(strcmp(ulogData.TopicNames,'yaw_estimator_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(yaw_estimator_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'yaw_estimator_status';
        disp('yaw_estimator_status')
    end

    % Make a table with the original sample frequency of each topic:
    disp('************************************************************')
    disp(' ')
    disp('Original topic sample rates:')
    fs_vec = fs_vec(2:end);
    sampleRatesTable = table(fs_vec', 'VariableNames', {'Rate [Hz]'}, 'RowNames', topicNames_vec);
    disp(sampleRatesTable)

    disp('************************************************************')
    
    % Find which uORB messages if any are not supported by this script
    unsupported_uORB = setdiff(dataTopicNames, topicNames_vec);
    if ismember('actuator_outputs',unsupported_uORB)
        idx = ismember(unsupported_uORB, 'actuator_outputs');
        unsupported_uORB(idx) = [];
    end
    if ~isempty(unsupported_uORB)
        % If 'uniqueNames' has data, display it
        disp(' ')
        disp('uORB topics not supported by this script:');
        disp(' ')
        disp(unsupported_uORB);    
        disp('************************************************************')
    end

    matFileName = ulog_name;
    matfile_path = strrep([ulog_path, ulog_name{1}], '.ulg', '.mat')
    
    % Create the directory if it doesn't exist
    if ~isfolder(fileparts(matfile_path))
        mkdir(fileparts(matfile_path));
    end

    ulog_parameters = struct();  % Create an empty structure for parameters
    output_vars = variableNames;
    output_vars{end+1} = 'PX4_parameters';
    save(matfile_path, output_vars{:});

    processTime = toc; % ends time

    fprintf('Processing completed in %.2f seconds.\n', processTime);
    fprintf('The .mat file path is:\n%s\n', matfile_path);    
end