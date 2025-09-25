clc;
close all;
clear variables;

%-------------------------------------------------------------------------%

% Reading Data
data = readstruct("processData.xml");

% Getting Data
T = data.T;
t = data.t;
Ks = data.Ks;

% Deviation Variables
Y = T - T(1);
X = t - t(1);

%-------------------------------------------------------------------------%

% Estimating tau_d

N = 50;
xdata = X(1:N);
ydata = Y(1:N);

p0 = [20, 1];
p = fmincon(@(p) errorFuncTauD(xdata, ydata, p), p0, ...
   [], [], [], [], [0, 0], [N, N]);

global tau_d
tau_d = p(1);
m = p(2);

Y_hat = zeros(size(xdata));

for i = 1:N
    if xdata(i) < p(1)
        Y_hat(i) = 0;
    else
        Y_hat(i) = p(2) * (xdata(i) - p(1));
    end
end

RMSE_tau_d = errorFuncTauD(xdata, ydata, p);
R2_tau_d = goodnessFit(ydata, Y_hat);

figure(1)
hold on
scatter(xdata, ydata, 50, "red", "filled", DisplayName="Recorded Data")
plot(xdata, Y_hat, LineWidth=1.75, Color="blue", DisplayName="Fitted Curve")
hold off
grid on
legend(Location="best")
xlabel("Time (sec)")
ylabel("Temperature (^oC)")
fontsize(20, "points")

%-------------------------------------------------------------------------%

% Estimating Kp
ydata = Y(end-N+1:end);
Y_inf = mean(ydata, "all");
Kp = Y_inf / Ks;

% Estimating tau_p
tau_p = Kp * Ks / m;

% Parameters
p0 = [Kp, tau_p];

% RMSE value
RMSE_unrefined = errorFunc(p0, Ks, X, Y);

% Estimated Output
Y_hat = firstOrderDelayTF(p0, Ks, X);

% Goodness of Fit
R2_unrefined = goodnessFit(Y, Y_hat);

%-------------------------------------------------------------------------%

% Refining Results

p = fmincon(@(p) errorFunc(p, Ks, X, Y), p0, ...
   [], [], [], [], [0, 0], []);

% Estimated Output
Y_hat = firstOrderDelayTF(p, Ks, X);

% Plotting figure
figure(2)
hold on
plot(X, Y, ...
    LineWidth=1.25, DisplayName="Recorded Data")
plot(X, Y_hat, ...
    LineWidth=1.75, DisplayName="Fitted Curve")
hold off
grid on
legend(Location="best")
axis([-50 950 -0.5 6])
xlabel("Time (sec)")
ylabel("Temperature (^oC)")
fontsize(20, "points")

% RMSE value
RMSE_refined = errorFunc(p, Ks, X, Y);

% Goodness of Fit
R2_refined = goodnessFit(Y, Y_hat);

%-------------------------------------------------------------------------%

% Writing Data to File
processModel = struct;
processModel.Kp = p(1);
processModel.tau_p = p(2);
processModel.tau_d = tau_d;

% writestruct(processModel, "processModel.xml")

%-------------------------------------------------------------------------%

% local functions

function [ Y_hat ] = firstOrderDelayTF(p, Ks, X)

% parameters
Kp = p(1);
tau_p = p(2);

global tau_d

% length of vector
n = length(X);
Y_hat = zeros(size(X));

% estimated value
for i = 1:n
    if X(i) < tau_d 
        Y_hat(i) = 0;
    else
        t_bar = X(i) - tau_d;
        Y_hat(i) = Kp .* Ks .* (1 - exp(-1 .* t_bar ./ tau_p));
    end
end 

end

function [ RMSE ] = errorFunc(p, Ks, X, Y)

% estimated value
Y_hat = firstOrderDelayTF(p, Ks, X); 

% error
err = Y - Y_hat;

% output
RMSE = sqrt(mean(err.^2, "all"));

end

function [ R2 ] = goodnessFit(Y, Y_hat)

Y_mean = mean(Y, "all");

SSR = sum((Y - Y_hat).^2, "all");
SST = sum((Y - Y_mean).^2, "all");

R2 = 1 - (SSR / SST);

end

function [ RMSE ] = errorFuncTauD(X, Y, p)

n = length(X);
Y_hat = zeros(size(X));

for i = 1:n
    if X(i) < p(1)
        Y_hat(i) = 0;
    else
        Y_hat(i) = p(2) * (X(i) - p(1));
    end
end

RMSE = sqrt(mean((Y - Y_hat).^2, "all"));

end

%-------------------------------------------------------------------------%