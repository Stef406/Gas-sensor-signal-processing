clc;
clear;
close all;
addpath('./data');

%% Import and preprocess data
raw = readtable('Steps_600 16-12-20.txt');
t_start = duration(inputdlg('Time of start operation'));                    % Prompt for start time
t_end = duration(raw.Var2(end));

% Remove data before start time and at end time
raw(raw.Var2 < t_start, :) = [];
raw(raw.Var2 == t_end, :) = [];
Conc = raw.Var6 * 1e-4;                                                     % Scale concentration
raw = [seconds(raw.Var2), Conc];                                            % Create raw time and concentration matrix

%% User input for valve operations
tup = seconds(duration(inputdlg('Time valve OPEN')));                       % Time valve opens
tdown = seconds(duration(inputdlg('Time valve CLOSED')));                   % Time valve closes
tel = 2 * 40; % Time extension
tf = tdown + tel; % Final time

% Find index values for up, down, and final times
idxu = find(raw(:, 1) == tup);
idxd = find(raw(:, 1) == tdown);
idxf = find(raw(:, 1) == tf);

% Extract raw data segments
Up = raw(idxu:idxd, 1);
Down = raw(idxd + 1:idxf, 1);
C_up = raw(idxu:idxd, 2);
C_down = raw(idxd + 1:idxf, 2);

% Normalize time relative to the start of "Up"
ti = Up - min(Up);
td = Down - min(Up);
s_up = [ti, C_up];
s_down = [td, C_down];

%% Filtering raw data
% Define window size for resampling (every 20 points)
N = numel(C_up);
idx = (0:20:N - 20).' + (1:20);
B = C_up(idx);
C_upfilt = mean(B, 2);                                                      % Filtered up data

M = numel(C_down);
idx = (0:20:M - 20).' + (1:20);
D = C_down(idx);
C_downfilt = mean(D, 2);                                                    % Filtered down data

% Combine time and filtered data
t = (0:length(C_upfilt) + length(C_downfilt) - 1)'; 
s_upfilt = [t(1:length(C_upfilt)), C_upfilt];
s_downfilt = [(t(length(C_upfilt) + 1):t(end))', C_downfilt];

%% Apply moving average and normalization
C_upavg = movmean(C_upfilt, 5);
C_upavg = C_upavg - min(C_upavg);                                           % Normalize

C_downavg = movmean(C_downfilt, 5);
C_downavg = C_downavg - min(C_downavg);                                     % Normalize

% Combine time and averaged data
s_upavg = [t(1:length(C_upavg)), C_upavg];
s_downavg = [t(1:length(C_downavg)), C_downavg];

%% Plot results
figure(1);
plot(s_up(:, 1), s_up(:, 2), s_down(:, 1), s_down(:, 2));                   % Raw data
hold on;
plot(s_upfilt(:, 1), s_upfilt(:, 2), s_downfilt(:, 1), s_downfilt(:, 2));   % Filtered data
plot(s_upavg(:, 1), s_upavg(:, 2), 'r-', 'LineWidth', 3);                   % Averaged up data
plot(s_downavg(:, 1), s_downavg(:, 2), 'k-', 'LineWidth', 3);               % Averaged down data
xlim([ti(1), td(end)]);
xlabel('Time (s)');
ylabel('CO_2 (%v/v)');
legend('Raw Up', 'Raw Down', 'Filtered Up', 'Filtered Down', 'Averaged Up', 'Averaged Down');
hold off;

%% Non-linear fit for up and down data
beta0 = [1 1 1 1];                                                          % Initial guess for parameters

% Fit for the "up" section using non-linear model fitting
mdlup = fitnlm(s_upavg(:, 1), s_upavg(:, 2), @cstr_up, beta0);
fit_up = cstr_up(mdlup.Coefficients.Estimate, s_upavg(:, 2));
tau_up = mdlup.Coefficients{4, 1};                                          % Response time in UP section of data

% Fit for the "down" section using non-linear model fitting
mdldown = fitnlm(s_downavg(:, 1), s_downavg(:, 2), @cstr_down, beta0);
fit_down = cstr_down(mdldown.Coefficients.Estimate, s_downavg(:, 2));
tau_down = mdldown.Coefficients{4, 1};                                      % Response time in DOWN section of data

tau = mean([tau_up, tau_down]);                                             % Average response time of the sensor

% Plot fitted results
figure(2);
plot(s_upavg(:, 1), s_upavg(:, 2), s_upavg(:, 1), fit_up);
hold on;
plot(s_downavg(:, 1), s_downavg(:, 2), s_downavg(:, 1), fit_down);
xlabel('Time (s)');
ylabel('CO_2 (%v/v)');
legend('Measured Up', 'Fitted Up', 'Measured Down', 'Fitted Down');
hold off;

%% Fitting functions for concentration profiles (system assumed as first order Continuous Stirred Tank Reactor, CSTR)
function fit_up = cstr_up(b, s_upavg)
    x = 0:length(s_upavg(:, 1)) - 1;
    fit_up = zeros(1, length(x));
    for k = 1:length(x)
        if x(k) <= b(3)
            fit_up(k) = b(1);
        else
            fit_up(k) = b(1) + b(2) * (1 - exp(-(x(k) - b(3)) / b(4)));
        end
    end
    fit_up = fit_up';
end

function fit_down = cstr_down(b, s_downavg)
    x = 0:length(s_downavg(:, 1)) - 1;
    fit_down = zeros(1, length(x));
    for k = 1:length(x)
        if x(k) <= b(3)
            fit_down(k) = b(1);
        else
            fit_down(k) = b(1) - b(2) * (1 - exp(-(x(k) - b(3)) / b(4)));
        end
    end
    fit_down = fit_down';
end
