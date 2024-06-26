clear; clc;

ulog_name = '20240327_Flight2.ulg';
ulog_path = 'C:\Users\gurt1\OneDrive - Virginia Tech\Documents\NSL Research\Flight Logs - HATCH\20240327\';
freq = 100; % [Hz]
cutoffFrequency = 10; % [Hz]

processUlogTimeSynced(ulog_path, ulog_name, freq, cutoffFrequency)