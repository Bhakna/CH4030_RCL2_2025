clc;
close all;
clear variables;

% Exp R1 : Biodiesel Synthesis
% RCL - 2

% Aayush Bhakna
% ch22b008

%-------------------------------------------------------------------------%

% Data from Experiment

t = [0, 15, 30, 45, 60, 75, 90, 105, 120];
AV = [4.15, 2.05, 0.56 3.59, 4.28, 4.49, 2.91, 4.62, 3.74];

%-------------------------------------------------------------------------%

% Parameter Estimation

p0 = [2, -1, -1];
% p0 = [0.1, 1, 1];
p = fmincon(@(p) errorFunc(AV, t, p), p0, [], []);

% Plotting Figure
t_data = linspace(0, 120, 100);

figure(1)
hold on
scatter(t, AV, 50, "red", "filled")
plot(t_data, AcidValue(t_data, p), "Color", "Blue", "LineWidth", 1.75)
hold off
grid on
legend('Data Points', 'Fitted Parameters', Location='best')
xlabel("Time (min)")
ylabel("Acid Value (mg KOH/g)")
fontsize(20, "points")

% Goodness of Fit
RMSE = errorFunc(AV, t, p);

AV_hat = AcidValue(t, p);
AV_mean = mean(AV, "all");

RSS = sum((AV - AV_hat).^2, "all");
TSS = sum((AV - AV_mean).^2, "all");

R2 = 1 - (RSS/TSS);

%-------------------------------------------------------------------------%

% Local function

function [ y ] = AcidValue(t, p)

m = p(1);
k = p(2);
c = p(3);

if m == 1
    y = c .* exp(k .* t);
else
    y = ((1-m) .* ((k.*t) + c)).^(1/(1-m));
end

end

function [ RMSE ] = errorFunc(AV, t, p)

AV_hat = AcidValue(t, p);
err = AV - AV_hat;

RMSE = sqrt(mean(err.^2, "all"));

end

%-------------------------------------------------------------------------%