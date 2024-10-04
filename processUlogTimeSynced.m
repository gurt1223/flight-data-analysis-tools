function processUlogTimeSynced(ulog_path, ulog_name, freq, CutoffFrequency)

% Convert ULOG to MAT file and Time Synchronize Data
%
% DESCRIPTION:
%   This script processes a ULog file by applying time synchronization and
%   filtering. It supports two modes of operation:
%   1. Interactive Mode: Without input arguments, it prompts the user to select
%      a ULOG file for conversion through a graphical interface.
%   2. Programmatic Mode: With provided ulog_path and ulog_name arguments, it directly
%      processes the specified ULOG file without user interaction.
%
% INPUTS:
%   ulog_path       - Path to the directory containing the ULog file (optional)
%   ulog_name       - Name of the ULog file (optional)
%   freq            - Desired frequency for time synchronization (Hz, optional)
%   cutoffFrequency - Cutoff frequency for the low-pass filter (Hz, optional)
%
% OUTPUTS:
%   A .mat file containing time-synchronized data
% 
% USAGE:
%   Interactive Mode: processUlogTimeSynced
%   Programmatic Mode: processUlogTimeSynced('C:/data', 'flightData.ulg')
%
% MATLAB REQUIREMENTS:
%   MATLAB R2023a or newer, UAV Toolbox
%
% WRITTEN BY:
%   Garrett D. Asper
%   Virginia Tech
%   Email: garrettasper@vt.edu
%
%   Patrick E. Corrigan
%   Virginia Tech
%   Email: patrickcorrigan@vt.edu
%
% HISTORY:
%   03 JUL 2024 - Created and debugged, GDA
%
% NOMENCLATURE:
%   The following topics will be time-synchronized if available in the log file:
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
% CREDIT:
%   Filtering elements of this code are thanks to Jeremy Hopwood from 
%   jwhgit01 FlightTest-Tools/src/FlightDataProcessing/FlightData.m

%% Input Validation and Error Handling

% ~~~~ Interactive Mode ~~~~ (Define these and then run the script)
    if nargin == 0 || isempty(ulog_name) || isempty(ulog_path) 
        clear; clc; % Use clear if you want to reset the environment in interactive mode
        % Prompt the user to select the Ulog file using the file explorer
        ulog_path = strrep(mfilename('fullpath'), mfilename, '');
        [ulog_name, ulog_path] = uigetfile([ulog_path, '*.ulg'], 'MultiSelect', 'on');
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

% Convert ulog_name to cell array if it's a single file selected
if ischar(ulog_name)
    ulog_name = {ulog_name};
end

% Validate the provided or selected file
if isequal(ulog_name,0) || isequal(ulog_path,0)
    error('File selection was cancelled or invalid file was provided.');
end

dt = 1 / freq; % Time step

%% Loop through each selected / inputted file
for fileIdx = 1:length(ulog_name)
tic % starts the timer

% Import Ulog
ulog_obj = ulogreader([ulog_path, ulog_name{fileIdx}]);

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

    %% actuator_controls_0
    if any(strcmp(ulogData.TopicNames,'actuator_controls_0'))
        fs_vec(end+1) = 1/median(diff(seconds(actuator_controls_0.timestamp)));
        topicNames_vec{end+1} = 'actuator_controls_0';
        actuator_controls_0_retimed = retime(actuator_controls_0, Time, 'pchip');
        actuator_controls_0 = actuator_controls_0_retimed;
    end

    %% actuator_outputs
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
    
    %% battery_status
    if any(strcmp(ulogData.TopicNames,'battery_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(battery_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'battery_status';

        % Retime
        battery_status_retimed = retime(battery_status, Time, 'nearest');
        battery_status = battery_status_retimed;
    end
    
    %% cpuload
    if any(strcmp(ulogData.TopicNames,'cpuload'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(cpuload.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'cpuload';

        % Retime
        cpuload_retimed = retime(cpuload, Time, 'pchip');
        cpuload = cpuload_retimed;
    end
    
    %% fcs_signals
    if any(strcmp(ulogData.TopicNames,'fcs_signals'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(fcs_signals.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'fcs_signals';

        % Retime
        fcs_signals = convertvars(fcs_signals,vartype('numeric'),'double');
        fcs_signals_retimed = retime(fcs_signals, Time, 'pchip');
        fcs_signals = fcs_signals_retimed;
    end

    %% fcs_signals_aux
    if any(strcmp(ulogData.TopicNames,'fcs_signals_aux'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(fcs_signals_aux.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'fcs_signals_aux';
    
        % Retime
        fcs_signals_aux = convertvars(fcs_signals_aux,vartype('numeric'),'double');
        fcs_signals_aux_retimed = retime(fcs_signals_aux, Time, 'pchip');
        fcs_signals_aux = fcs_signals_aux_retimed;
    end
   
      %% esc_status
    if any(strcmp(ulogData.TopicNames,'esc_status'))
        % Find the index of 'esc_status' in the TopicNames array
        idx = find(strcmp(ulogData.TopicNames, 'esc_status'));
        fs_vec(end+1) = 1/median(diff(seconds(esc_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'esc_status';
    
        % Access the corresponding TopicMessages using the index
        esc_status = ulogData.TopicMessages{idx};
        timestamp_ = seconds(seconds(esc_status.timestamp));
        
        % Convert all variables in esc_status to double
        allVars = vartype('numeric'); % Selects all numeric variables
        esc_status = convertvars(esc_status, allVars, 'double');
        
        % Retime the ESC signals
        esc_status_retimed = retime(esc_status, Time, 'pchip');
        esc_status = esc_status_retimed;
        

    end


    %% estimator_states
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
    
    %% estimator_status
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
    
    %% estimator_event_flags
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

    %% event
    if any(strcmp(ulogData.TopicNames,'event'))        
        % Retiming has not yet been incorporated for this topic
        disp('event')
    end
    
    %% input_rc
    if any(strcmp(ulogData.TopicNames,'input_rc'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(input_rc.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'input_rc';
    
        % Retime
        input_rc_retimed = retime(input_rc, Time, 'nearest');
        input_rc = input_rc_retimed;
    end
        
    %% sensor_combined
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

    
    %% sensor_gps
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

    %% vehicle_acceleration
    if any(strcmp(ulogData.TopicNames,'vehicle_acceleration'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_acceleration.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_acceleration';

        % Retime
        vehicle_acceleration_retimed = retime(vehicle_acceleration, Time, 'pchip');
        vehicle_acceleration = vehicle_acceleration_retimed;
    end
    
    %% vehicle_air_data
    if any(strcmp(ulogData.TopicNames,'vehicle_air_data'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_air_data.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_air_data';

        disp('vehicle_air_data')
    end
    
    %% vehicle_angular_velocity
    if any(strcmp(ulogData.TopicNames,'vehicle_angular_velocity'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_angular_velocity.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_angular_velocity';

        % Retime
        vehicle_angular_velocity_retimed = retime(vehicle_angular_velocity, Time, 'pchip');
        vehicle_angular_velocity = vehicle_angular_velocity_retimed;
    end
    
    %% vehicle_attitude
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
    
    %% vehicle_attitude_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_attitude_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_attitude_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_attitude_setpoint';

        % Retime
        vehicle_attitude_setpoint_retimed = retime(vehicle_attitude_setpoint, Time, 'previous');
        vehicle_attitude_setpoint = vehicle_attitude_setpoint_retimed;
    end
    
    %% vehicle_control_mode
    if any(strcmp(ulogData.TopicNames,'vehicle_control_mode'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_control_mode.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_control_mode';

        % Retime
        vehicle_control_mode_retimed = retime(vehicle_control_mode, Time, 'previous');
        vehicle_control_mode = vehicle_control_mode_retimed;
    end
    
    %% vehicle_global_position
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
    
    %% vehicle_local_position
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
    
    %% vehicle_local_position_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_local_position_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_local_position_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_local_position_setpoint';
        
        % Retime
        vehicle_local_position_setpoint_retimed = retime(vehicle_local_position_setpoint(:,["x","y","z","vx","vy","vz"]),Time,'previous');
        vehicle_local_position_setpoint = vehicle_local_position_setpoint_retimed;
    end
    
    %% vehicle_rates_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_rates_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_attitude.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_rates_setpoint';

        disp('vehicle_rates_setpoint')
    end
    
    %% vehicle_status
    if any(strcmp(ulogData.TopicNames,'vehicle_status'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_status.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_status';
        
        % Retime
        vehicle_status_retimed = retime(vehicle_status, Time, 'nearest');
        vehicle_status = vehicle_status_retimed;
    end
    
    %% vehicle_torque_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_torque_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_torque_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_torque_setpoint';
        disp('vehicle_torque_setpoint')
    end
    
    %% vehicle_thrust_setpoint
    if any(strcmp(ulogData.TopicNames,'vehicle_thrust_setpoint'))
        % Find the original sample rate
        fs_vec(end+1) = 1/median(diff(seconds(vehicle_thrust_setpoint.timestamp))); % Sampling frequency
        topicNames_vec{end+1} = 'vehicle_thrust_setpoint';
        disp('vehicle_thrust_setpoint')
    end
    
    %% yaw_estimator_status
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

    matFileName = ulog_name{fileIdx};
    matfile_path = strrep([ulog_path, ulog_name{fileIdx}], '.ulg', '.mat');
    
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
end