clc;
clear variables;
close all;

%-------------------------------------------------------------------------%

% Data Manipulation

givendata = readtable("RCL2_Aug06_2025.xlsx");
givendata = givendata(:, ["Time", "Untitled", "Untitled1", "Untitled2"]);
givendata = table2cell(givendata);

data = zeros(length(givendata), 4);
data(:, 2:4) = cell2mat(givendata(:, 2:4));

for i = 1:length(givendata)
    startTime = givendata(i, 1);
    startTime = char(string(startTime));
    startTime = floor(seconds(timeofday(datetime(startTime(13:end)))));
    data(i, 1) = startTime;
end

data = data(2:end, :);
data(:, 1) = data(:, 1) - data(1, 1);

clear givendata i startTime;
clc;

% Data Information :
%   Column 1 : "Time in sec"
%   Column 2 : "Set Point in ^oC"
%   Column 3 : "Temperature in ^oC"
%   Column 4 : "Volume Flow Rate in %"

%-------------------------------------------------------------------------%

% % Set Point : 70
% setPoint = 70;
% indexStart = 4040;
% indexEnd = 4122;
% 
% % Normalized Data
% vec_time = data(indexStart:indexEnd, 1) - data(indexStart, 1);
% vec_temp = data(indexStart:indexEnd, 3) - data(indexStart, 3);
% setPoint = setPoint - data(indexStart, 3);
% 
% % Process Function
% p0 = [10, 10, 1.5];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [0, 0, 0];
% ub = [];
% params_70 = fmincon(@(pval) minFunc(pval, vec_time, vec_temp, setPoint), ...
%     p0, A, b, Aeq, beq, lb, ub);
% 
% % Calculated Value
% s = tf('s');
% G = @(p) p(1) / ( ((p(2) * s)^2) + (2 * p(2) * p(3) * s) + 1 );
% 
% temp_hat = setPoint .* step(G(params_70), 0:vec_time(end));
% 
% % Result Plot
% figure(2)
% hold on
% plot([vec_time(1), vec_time(end)], setPoint .* [1, 1], ...
%     LineWidth=2, Color="Blue")
% scatter(vec_time, vec_temp, 30, "red", "filled")
% plot(0:vec_time(end), temp_hat, ...
%     LineWidth=2, Color="Black")
% hold off
% grid on
% axis([0 400 0 6])
% legend('Set Point', 'Observed Data', 'Fitted Data', Location='best')
% xlabel('Normalized Time (sec)')
% ylabel('Normalized Temperature (celcius)')
% title('Set Point : 70 Celcius')
% fontsize(18, "points")

%-------------------------------------------------------------------------%

% Set Point : 75
setPoint = 75;
indexStart = 4676;
indexEnd = 4742;

% Normalized Data
vec_time = data(indexStart:indexEnd, 1) - data(indexStart, 1);
vec_temp = data(indexStart:indexEnd, 3) - data(indexStart, 3);
setPoint = setPoint - data(indexStart, 3);

% Process Function
p0 = [10, 10, 1.5];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0, 0, 0];
ub = [];
params_75 = fmincon(@(pval) minFunc(pval, vec_time, vec_temp, setPoint), ...
    p0, A, b, Aeq, beq, lb, ub);

% Calculated Value
s = tf('s');
G = @(p) p(1) / ( ((p(2) * s)^2) + (2 * p(2) * p(3) * s) + 1 );

temp_hat = setPoint .* step(G(params_75), 0:vec_time(end));

% Result Plot
figure(2)
hold on
plot([vec_time(1), vec_time(end)], setPoint .* [1, 1], ...
    LineWidth=2, Color="Blue")
scatter(vec_time, vec_temp, 30, "red", "filled")
plot(0:vec_time(end), temp_hat, ...
    LineWidth=2, Color="Black")
hold off
grid on
% axis([0 400 0 6])
legend('Set Point', 'Observed Data', 'Fitted Data', Location='best')
xlabel('Normalized Time (sec)')
ylabel('Normalized Temperature (celcius)')
title('Set Point : 75 Celcius')
fontsize(18, "points")

%-------------------------------------------------------------------------%

% min function

function [ output ] = minFunc(pval, vec_time, vec_temp, setPoint)

s = tf('s');
G = @(p) p(1) / ( ((p(2) * s)^2) + (2 * p(2) * p(3) * s) + 1 );

temp_hat = setPoint .* step(G(pval), 0:vec_time(end));

vec_temp_hat = temp_hat(vec_time + 1);
output = sum((vec_temp_hat - vec_temp).^2, "all");

end

%-------------------------------------------------------------------------%