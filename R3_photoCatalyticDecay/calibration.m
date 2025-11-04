clc;
close all;
clear variables;

%-------------------------------------------------------------------------%

% Recorded Data

% concentration (in ppm)
C = [0, 5, 10, 15, 20, 25];

% absorbance
A = [0, 0.386422, 0.759161, 1.12572, 1.48722, 1.90242];

% max absorbing wavelength (in nm)
lambda_max = 463.5; 

%-------------------------------------------------------------------------%

% Curve Fitting

% slope
C_mean = mean(C, "all");
A_mean = mean(A, "all");
m = sum((C - C_mean).*(A - A_mean), "all") / sum((A - A_mean).^2, "all");

% polynomial form (no intercept)
p = [m, 0];

% plot
figure(1)
hold on
scatter(A, C, ...
    75, "red", "filled", ...
    DisplayName="Data Points")
plot(A, polyval(p, A), ...
    "Color", "blue", "LineWidth", 1.5, "Marker", "x", ...
    DisplayName="Fitted Line")
hold off
grid minor
legend(Location="best")
xlabel('Absorbance at 463.5 nm')
ylabel('Concentration (ppm)')
fontsize(20, "points")

% results
C_hat = polyval(p, A);
RSS = sum((C - C_hat).^2, "all");
TSS = sum((C - C_mean).^2, "all");
RESULTS_R2 = 1 - (RSS / TSS);
RESULTS_RMSE = sqrt(mean((C - C_hat).^2, "all"));

%-------------------------------------------------------------------------%