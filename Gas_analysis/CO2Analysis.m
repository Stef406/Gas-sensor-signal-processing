clc; clear; close all;
format short
addpath('./data');

%% Import gas data and pre-process
% Get start time of operation from user input
t_start = duration(inputdlg('Time of start operation'));

% Scaling factor for CO2 sensor data (depends on sensor range)
s_CO2 = 10;

% Read raw CO2 data from text file
raw_CO2 = readtable('Gas 8-12-20_600.txt', 'PreserveVariableNames', true);

% Process the CO2 data
raw_CO2 = dioxide(raw_CO2, s_CO2, t_start);

%% Postprocessing parameters
f_diox = 20;                                                                % Sampling frequency (Hz) for CO2 sensor
n_diox = 20;                                                                % Window width for moving average
Q = 1;                                                                      % Sampling volume flow rate (l/min)
T = 25 + 273.15;                                                            % Temperature in sensor (K)
P = 101325 * 10^-4;                                                         % Pressure in sensor (atm)

% Get CO2 inlet times
inlet = duration(inputdlg({'time_0', 'elapsed'}, 'CO2 in', [1 50; 1 50]));

% Convert CO2 inlet times times to seconds
time = [seconds(inlet(1)), seconds(inlet(2))];
tin = time(1,1) + time(1,2);

% Perform deconvolution on the CO2 data
[C_mdiox, C_tdiox, mol_mdiox, mol_tdiox] = deconv(raw_CO2, tin, f_diox, n_diox, Q, T, P);

%% Plot results
figure(1)
plot(C_mdiox(:,1), C_mdiox(:,2), C_tdiox(:,1), C_tdiox(:,2));
xlabel('Time (s)');
ylabel('CO_2 (%v/v)');
yline(0); % Add y-axis line at 0
legend('Measured', 'Deconvoluted');

%% Preprocessing function for CO2 data
function Conc_diox = dioxide(data, scale, t_start)
    % Set variable names
    data.Properties.VariableNames = {'date', 'time', 'Z', 'Cfilt', 'z', 'Cunfilt'};
    
    % Get end time of CO2 data
    t_end_CO2 = duration(data.time(end));
    
    % Filter data between start and end times
    data(data.time < t_start, :) = [];
    data(data.time == t_end_CO2, :) = [];
    
    % Scale CO2 concentration (volumetric fraction %)
    CO2 = data.Cfilt * scale * 1e-4;
    
    % Return time in seconds and scaled CO2 concentration
    Conc_diox = [seconds(data.time), CO2];
end

%% Deconvolution function
function [C_m, C_t, mol_m, mol_t] = deconv(data, tin, freq, avg, Q, Temp, Press)
    % Set duration of the concentration profile of interest
    tel = 2 * 230;                                                          % Total time elapsed (s)
    tf = tin + tel;
    
    % Find the index corresponding to start and end times
    idxu = find(data(:, 1) == tin);
    idxf = find(data(:, 1) == tf);
    
    % Extract raw CO2 data within the time range of interest
    C_raw = data((idxu:idxf), 2);
    
    % Resample data based on the specified sampling frequency
    if freq > 1
        N = numel(C_raw);
        idx = (0:freq:N - freq).' + (1:freq);
        B = C_raw(idx);
        C_filt = mean(B, 2);                                                % Apply resampling by averaging
    else
        C_filt = C_raw;
    end
    
    % Generate time vector
    t = (0:length(C_filt) - 1)';
    
    % Apply moving average filter to the data
    if ~isempty(avg)
        C_avg = movmean(C_filt, avg);
    else
        C_avg = C_filt;
    end
    
    % Normalize concentration data by subtracting the minimum value
    C_avg = C_avg - min(C_avg);
    
    % Measured concentration profile (post-processed)
    C_m = [t(1:length(C_avg)), C_avg];
    
    % Deconvolution: Estimate the sensor response time and adjust data
    tau = 13.70;                                                            % CO2 sensor response time (depends on the experimental conditions, determined from StepDeconv.m)
    
    % Apply deconvolution to find the true concentration profile (assumed
    % first order system)
    C_t = C_m(2:end, 2) + tau .* diff(C_m(:, 2)) ./ diff(C_m(:, 1));
    C_t = [t(1:length(C_t)), C_t];
    C_t(:, 2) = smooth(C_t(:, 2));                                          % Deconvoluted data
    
    %% Calculate moles of gas using the ideal gas law
    R = 0.0821;                                                             % Universal gas constant (l atm / mol K)
    
    % Calculate moles for measured and deconvoluted data
    mol_m = (Press * Q / (R * Temp)) * trapz(C_m(:, 1), C_m(:, 2)) / (100 * 60);
    mol_t = (Press * Q / (R * Temp)) * trapz(C_t(:, 1), C_t(:, 2)) / (100 * 60);
end
