% % Define necessary variables
% ulog_name = '20240324_Flight1.ulg';
% ulog_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240324\';
% freq = 60; % [Hz] (also results in frames per second of vide)
% cutoffFrequency = 10; % for filtering angular rates
% trim_start = 0; % start time of the video corresponding to flight log time (seconds)
% trim_end = 0; % end time of the video corresponding to flight log time (seconds)
%               % leave 0 for full flight log
% mp4FileName = '20240324_Flight1'; % name of the visualization video file (do not append with .mp4)
% 
% % Call the function
% animateFlight(ulog_name, ulog_path,freq, cutoffFrequency, trim_start, trim_end, mp4FileName)

% Define necessary variables
ulog_name = '20240324_Flight2.ulg';
ulog_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240324\';
freq = 60; % [Hz] (also results in frames per second of vide)
cutoffFrequency = 10; % for filtering angular rates
trim_start = 50; % start time of the video corresponding to flight log time (seconds)
trim_end = 70; % end time of the video corresponding to flight log time (seconds)
              % leave 0 for full flight log
mp4FileName = '20240324_Flight2_NumMethFinalProject'; % name of the visualization video file (do not append with .mp4)

% Call the function
animateFlight(ulog_name, ulog_path,freq, cutoffFrequency, trim_start, trim_end, mp4FileName)