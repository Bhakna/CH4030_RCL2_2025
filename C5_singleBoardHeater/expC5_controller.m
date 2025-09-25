clc;
clear variables;
close all;

%-------------------------------------------------------------------------%

% Get parameters (from "expC5_firstOrder_delay.m")
data = readstruct("processModel.xml");

% Model estimates
Kp = data.Kp;
tau_p = data.tau_p;
tau_d = data.tau_d;

% No PID
p0 = [1, 1e12, 0];

%-------------------------------------------------------------------------%

% Closed Loop

% Calculating ultimate values
x0 = [100, 0.1];  % guess
options = optimoptions("fmincon", ...
    "StepTolerance", 1e-6, ...
    "ConstraintTolerance", 1e-9, ...
    "Display", "none");
x = fmincon(@(x) minFunc(Kp, tau_p, tau_d, x), x0, ...
    [], [], [], [], [1e-6, 1e-6], [], [], options);

% Ultimate Gain
Kcu = x(1);

% Ultimate Frequency / Period
wcu = x(2);
Pcu = 2 * pi / wcu;

% Closed-Loop ZN tuning
Kc = Kcu / 1.7;
tau_I = Pcu / 2;
tau_D = Pcu / 8;

% Plotting Tuned Process
p = [Kc, tau_I, tau_D];
T = linspace(0, 1000, 1000); % sec
Y = step(systemTF(Kp, tau_p, tau_d, p0), T);
Y_tuned = step(systemTF(Kp, tau_p, tau_d, p), T);
figure(1)
hold on
plot(T, Y, Color="blue", LineWidth=1.5, DisplayName="No Tuning")
plot(T, Y_tuned, Color="red", LineWidth=1.5, DisplayName="ZN Tuning")
hold off
grid minor
legend(Location="best")
xlabel("Time (sec)")
ylabel("Temperature (^oC)")
fontsize(20, "points")

% Writing to struct
results = struct;
results.Closed = p;

%-------------------------------------------------------------------------%

% Open Loop

% Open-Loop ZN tuning
Kc = 1 / Kp;
tau_I = 1e12;
tau_D = 0;

% Plotting Tuned Process
p = [Kc, tau_I, tau_D];
T = linspace(0, 1000, 1000); % sec
Y = step(systemTF(Kp, tau_p, tau_d, p0, "open"), T);
Y_tuned = step(systemTF(Kp, tau_p, tau_d, p, "open"), T);
figure(2)
hold on
plot(T, Y, Color="blue", LineWidth=1.5, DisplayName="No Tuning")
plot(T, Y_tuned, Color="red", LineWidth=1.5, DisplayName="ZN Tuning")
hold off
grid minor
legend(Location="best")
xlabel("Time (sec)")
ylabel("Temperature (^oC)")
fontsize(20, "points")

% Writing to struct
results.Open = p;

%-------------------------------------------------------------------------%

% process transfer function
function [ G_p ] = processTF(Kp, tau, tau_d)

% process 
s = tf('s');
G_p = Kp * exp(-1 * tau_d * s) / (1 + (tau * s));

end

% overall transfer function
function [ G_overall ] = systemTF(Kp, tau, tau_d, p, type)

% default value
if nargin < 5
    type = "closed";
end

% process 
s = tf('s');
G = processTF(Kp, tau, tau_d);

% controller parameters
Kc = p(1);
tau_I = p(2);
tau_D = p(3);

% controller
C = Kc * (1 + (1/(tau_I * s)) + (tau_D * s));

% output
if type == "closed" 
    G_overall = G * C / (1 + (G * C));
else
    G_overall = G * C;
end

end

% function to calculate Ultimates
function [ RMSE ] = minFunc(Kp, tau_p, tau_d, x)

% input
Kcu = x(1);
wcu = x(2);

% function
f_x = zeros(2, 1);
f_x(1) = 1 + (Kp * Kcu * cos(tau_d * wcu));
f_x(2) = (tau_p * wcu) - (Kp * Kcu * sin(tau_d * wcu));

% output
RMSE = sqrt(mean(f_x.^2, "all"));

end

%-------------------------------------------------------------------------%