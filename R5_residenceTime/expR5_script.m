clc;
clear variables;
close all;

% CH4030 - Reaction and Control Lab
% Group - 2

% Aayush Bhakna
% CH22B008

%-------------------------------------------------------------------------%

% reading recorded data

data = readstruct("ExpR5_dataRecords.xml");
currentCase = data.Case1;

% Time
t = currentCase.T;
delT = 0.5;

% Delta Time
t_data = linspace(min(t, [], "all"), max(t, [], "all"), 100);

% Conductance
G = currentCase.C;
G_0 = G(1);
G_bar = G - G_0;

% clearing error
for i = 1:length(G_bar)
    if G_bar(i) < 0
        G_bar(i) = 0;
    end
end

% Residence Time Distribution
E_t = G_bar ./ (sum(G_bar, "all") * delT);

param0 = [4, 1];
param = fmincon(@(p) errorFunc(t, E_t, p), param0, [], []);

figure(1)
hold on
plot(t_data, normalDist(t_data, param), Color="blue", LineWidth=1.75)
scatter(t, E_t, 50, "red", "filled")
hold off
grid on
legend('Fitted Normal Dist', 'Collected Data', Location='northeast')
xlabel("t (min)")
ylabel("E(t)")
title("Residence Time Distribution")
fontsize(20, "points")

%-------------------------------------------------------------------------%

% Calculating Peclet Number

% Mean Residence Time
tm = sum((t .* E_t .* delT), "all");

% Variance
var = sum((((t - tm).^2) .* E_t .* delT), "all");

% Peclet Number
Pe = fsolve(@(x) fsolveFunc(x, tm, var), 1);

% Creating result table
results = table;
results.Time = t';
results.G_bar = G_bar';
results.E_t = E_t';
results.t_G_bar = (t .* G_bar)';
results.t2_G_bar = ((t.^2) .* G_bar)';

disp(results)

%-------------------------------------------------------------------------%

% reading recorded data

function [ y ] = normalDist(x, p)

mu = p(1);
var = p(2);
y = exp(-1 .* ((x - mu).^2) ./ (2 * var)) ./ sqrt(2 * pi * var);

end

function [ RMSE ] = errorFunc(x, y, p)

y_hat = normalDist(x, p);
err = y - y_hat;

RMSE = sqrt(mean(err.^2, "all"));

end

function [ f_x ] = fsolveFunc(x, tm, var)

f0 = var / (2 * (tm^2));
f_x = ((x - 1 + exp(-x)) / (x^2)) - f0;

end

%-------------------------------------------------------------------------%