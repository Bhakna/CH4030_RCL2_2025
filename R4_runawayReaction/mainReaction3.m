clc;
close all;
clear variables;

data = readmatrix("expR4_recordedData.xlsx");
t = data(:, 1);
T = data(:, 5);

global UA; global T_inf; global V_mix; 
global MD; global CvD;

UA = 0.911; % J/kg C
T_inf = 40.4; % celcius

V_mix = 100 * 1e-6; % m^3

global CA0; global CB0; global CC0; global CD0;
global Tin;

Tin = T(1);

MA = 0.248; % kg/mol
rhoA = 1667; % kg/m^3
CA0 = 0.45; % mol/m^3

MB = 0.034;
rhoB = 1450;
CB0 = 0.90;

MC = 0.142;
rhoC = 2664;
CC0 = 0;

MD = 0.018;
rhoD = 1000;
alpha = 1 - (CA0 * MA / rhoA) - (CB0 * MB / rhoB);
CD0 = rhoD * alpha / MD;
CvD = 4180; % J/kg*K

P_guess = [-40, -1e5];
P_calc = fmincon(@(P) errorFunc(P, t, T), P_guess, ...
    [], []);

P = P_calc;
lnA = P(1);
delHr = P(2);
X0 = [CA0; CB0; CC0; CD0; Tin];
[~, X] = ode23s(@(t, X) odeFunc(t, X, lnA, delHr), t, X0);
CA = X(:, 1);
CB = X(:, 2);
CC = X(:, 3);
CD = X(:, 4);
T_hat = X(:, 5);

figure(1)
hold on
plot(t, CA, LineWidth=1.5, LineStyle="-", DisplayName="Thiosulphate")
plot(t, CB, LineWidth=1.5, LineStyle="-", DisplayName="Hydrogen Peroxide")
plot(t, CC, LineWidth=1.5, LineStyle="-", DisplayName="Sulphate")
hold off
grid on
legend(Location="best")
xlabel("Time (sec)")
ylabel("Concentration (mol/m^3)")
fontsize(20, "points")

figure(2)
hold on
plot(t, T, LineWidth=1.5, LineStyle="-", DisplayName="Data Points")
plot(t, T_hat, LineWidth=1.5, LineStyle="-.", DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel("Time (sec)")
ylabel("Temperature (celcius)")
fontsize(20, "points")

RMSE = sqrt(mean((T - T_hat).^2, "all"));
T_mean = mean(T, "all");
R2 = 1 - (( sum((T - T_hat).^2, "all") ) / ( sum((T - T_mean).^2, "all") ));

function [ k ] = rateConst(A, del_Hr, T)

T_K = T + 273.15;

R = 8.314; % J/mol*K
k = exp( A - (del_Hr / (R * T_K)) );

end

function [ f ] = odeFunc(t, X, A, del_Hr)

% global variables
global UA; global T_inf; global V_mix; 
global MD; global CvD;

% input
CA = X(1);
CB = X(2);
CC = X(3);
CD = X(4);
T = X(5);

% rate
r = rateConst(A, del_Hr, T) * CA * CB;

% coefficients
delT = T - T_inf;
rhoCvV = V_mix * CD * MD * CvD;

% output
f = zeros(5, 1);
f(1) = -1 * r;
f(2) = -4 * r;
f(3) = 2 * r;
f(4) = 3 * r;
f(5) = -1 * ( (del_Hr * r) + (UA * delT) ) / rhoCvV;

end

function [ RMSE ] = errorFunc(P, t, T)

A = P(1);
del_Hr = P(2);

global CA0; global CB0; global CC0; global CD0;
global Tin;

X0 = [CA0; CB0; CC0; CD0; Tin];
[~, X] = ode23s(@(t, X) odeFunc(t, X, A, del_Hr), t, X0);

T_hat = X(:, 5);
err = T - T_hat;

RMSE = sqrt(mean(err.^2, "all"));

end
