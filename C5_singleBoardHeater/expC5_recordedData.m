clc;
close all;
clear variables;

%-------------------------------------------------------------------------%

% Reading Data
data = readmatrix("rcl2.txt");

% Getting Data
serialNum = data(:, 1);    % index (-)
CO = data(:, 2);           % controller output (%)
D = data(:, 3);            % disturbance (%)
PV = data(:, 4);           % process variables (^oC)
timeStamp = data(:, 5);    % instantaneous CPU time (ms)

% Converting Units of Data
timeStamp = (timeStamp - timeStamp(1)) ./ 1000; % seconds

%-------------------------------------------------------------------------%

% Plotting figure
figure(1)

subplot(3, 1, 1)
plot(timeStamp, CO, LineWidth=1.25, Color="blue")
grid on
xlabel("Time (sec)")
ylabel("CO (%)")

subplot(3, 1, 2)
plot(timeStamp, D, LineWidth=1.25, Color="black")
grid on
xlabel("Time (sec)")
ylabel("D (%)")

subplot(3, 1, 3)
plot(timeStamp, PV, LineWidth=1.25, Color="red")
grid on
xlabel("Time (sec)")
ylabel("PV (^oC)")

fontsize(20, "points")

%-------------------------------------------------------------------------%

% Data after Step-Change
processData = struct;
processData.t = timeStamp(101:end);
processData.T = PV(101:end);
processData.Ks = CO(101) - CO(100);

% Write Data to File
% writestruct(processData, "processData.xml")

%-------------------------------------------------------------------------%