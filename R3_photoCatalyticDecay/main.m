clc;
close all;
clear variables;

%-------------------------------------------------------------------------%

% Recorded Data

% time (in min)
t = 0:10:90;

% absorbance at max wavelength
A_max = [1.75159, 1.49905, 1.27180, 1.09803, 0.984364, ...
    0.847860, 1.01707, 0.953266, 0.987498, 0.749132];

% max absorbing wavelength (in nm)
lambda_max = [461.2, 457.6, 452.8, 445.0, 403.0, ...
    402.2, 341.0, 342.8, 340.4, 336.4]; 

% absorbance at fixed wavelength
A = [1.75065, 1.49022, 1.25039, 1.06117, 0.920248, ...
    0.735692, 0.655930, 0.574247, 0.456763, 0.355654];

% fixed absorbing wavelength (in nm)
lambda = 463.0; 

% slope of variable conversion (A -> C)
m = 13.2720;

% plot
figure(1)
hold on
yyaxis left
plot(t, lambda_max, ...
    "Color", "blue", "LineWidth", 1.5, "Marker", "o", ...
    DisplayName="Wavelength")
ylabel('Max Absorbing Wavelength (nm)')
yyaxis right
plot(t, A_max, ...
    "Color", "red", "LineWidth", 1.5, "Marker", "*", ...
    DisplayName="Absorbance")
hold off
grid minor
legend(Location="best")
xlabel('Time (min)')
ylabel('Max Absorbance')
fontsize(20, "points")

% plot
figure(2)
hold on
plot(lambda_max, A_max, ...
    "Color", "red", "LineWidth", 1.5, "Marker", "o", ...
    DisplayName="Data Points")
hold off
grid minor
legend(Location="best")
xlabel('Max Absorbing Wavelength (nm)')
ylabel('Max Absorbance')
fontsize(20, "points")

%-------------------------------------------------------------------------%

% Calculating Conc

% concentration at fixed wavelength
C = A .* m;

% plot
figure(3)
hold on
plot(t, C, ...
    "Color", "blue", "LineWidth", 1.5, "Marker", "o", ...
    DisplayName="Data Points")
hold off
grid minor
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Concentration at 463.0 nm (ppm)')
fontsize(20, "points")

%-------------------------------------------------------------------------%

% DIFFERENTIAL ANALYSIS

% derivative of concentration in (ppm/min)
dCdt = zeros(size(C));

% forward diff
dCdt(1) = (C(2) - C(1)) / (t(2) - t(1));

% central diff
for i = 2:length(t)-1
    dCdt(i) = (C(i+1) - C(i-1)) / (t(i+1) - t(i-1));
end

% backward diff
dCdt(end) = (C(end) - C(end-1)) / (t(end) - t(end-1));

% fitting polynomial
xdata = log10(C);
ydata = log10(-dCdt);
D_p = polyfit(xdata, ydata, 1);

% plotting derivative vs conc in log-log
% x = linspace(-0.6, 0.6, 100);
figure(4)
hold on
scatter(xdata, ydata, 50, "red", "filled", DisplayName="Data Points")
plot(xdata, polyval(D_p, xdata), "Color", "blue", "LineWidth", 1.75, DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel('Concentration (ppm)')
ylabel('-1 * Conc Derivative (ppm/min)')
fontsize(20, "points")

% Parameter values
D_a = D_p(1);
D_kr = 10^(D_p(2));

% Goodness of fit
D_r2 = goodnessFit(ydata, polyval(D_p, xdata));

%-------------------------------------------------------------------------%

% INTEGRAL ANALYSIS

INT_C = zeros(size(C));

for i = 2:length(t)
    INT_C(2:i) = INT_C(2:i) + ((C(i) .^ (-D_a)) * (C(i) - C(i-1)));
end

INT_C = INT_C(2:end);
INT_p = polyfit(t(2:end), INT_C, 1);

% plotting integral vs time
% x = linspace(0, 55, 100);
figure(5)
hold on
scatter(t(2:end), INT_C, 50, "red", "filled", DisplayName="Data Points")
plot(t(2:end), polyval(INT_p, t(2:end)), "Color", "blue", "LineWidth", 1.75, DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Integral')
fontsize(20, "points")

INT_kr = INT_p(1);

% Goodness of fit
INT_r2 = goodnessFit(INT_C, polyval(INT_p, t(2:end)));

%-------------------------------------------------------------------------%

% REGRESSION ANALYSIS

% initial concentration (in ppm)
C0 = C(1);

% parameters from differential analysis (taken as guess here)
p0 = [D_kr, D_a];

% estimating (fitting) values using regression
options = optimoptions("lsqnonlin", ...
    "Algorithm", "levenberg-marquardt", ...
    "FunctionTolerance", 1e-8, ...
    "MaxIterations", 1e4, ...
    "Display", "none");
p = lsqnonlin(@(p) errorFunc(C0, C, t, p), ...
    p0, [], [], [], [], [], [], [], options);

% solving ode
[~, C_REG] = ode45(@(t, C) reaxEq(C, t, p), t, C0);

% plotting graph
figure(6)
hold on
scatter(t, C, 50, "red", "filled", DisplayName="Data Points")
plot(t, C_REG, Color="blue", LineWidth=1.75, Marker="o", DisplayName="REG Fitting")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Concentration (ppm)')
fontsize(20, "points")

% Final Error
REG_RMSE = errorFunc(C0, C, t, p);

% Goodness of fit
REG_r2 = goodnessFit(C, (C_REG)');


%-------------------------------------------------------------------------%

% LOCAL FUNCTIONS

% Model Equations (ODEs)
function [ dCdt ] = reaxEq(C, t, p)

% parameters
kr = -p(1);
a = p(2);

% output
dCdt = kr * (C ^ a);

end

% Minimizing Function 
function [ RMSE ] = errorFunc(C0, C, t, p)

% solving ode
[~, ode_C] = ode45(@(t, C) reaxEq(C, t, p), t, C0);
C_hat = ode_C;

% error
err = C' - C_hat;

% output
RMSE = sqrt(mean(err.^2, "all"));

end

% Goodness of Fit
function [ R2 ] = goodnessFit(Y, Y_hat)

Y_mean = mean(Y, "all");
RSS = sum((Y - Y_hat).^2, "all");
TSS = sum((Y - Y_mean).^2, "all");

R2 = 1 - (RSS / TSS);

end

%-------------------------------------------------------------------------%

% LICENSE IN NAME : Aayush Bhakna (ch22b008)

%-------------------------------------------------------------------------%