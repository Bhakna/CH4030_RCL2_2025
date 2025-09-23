clc;
close all;
clear variables;

% Experiment - R2
% PHASE TRANSFER CATALYSIS

% Aayush Bhakna
% ch22b008

%-------------------------------------------------------------------------%

% NOTATIONS USED

% NC -> no catalyst
% C -> with catalyst

% A -> NaOH
% B -> CH3COO-C4H9

% 1 -> Aqueous Phase
% 2 -> Organic Phase

%-------------------------------------------------------------------------%

% INITIALIZING DATA

% Recorded Data
t = 0:5:50;
delV_NC = [6.0, 3.7, 5.4, 5.1, 4.8, 4.1, 5.3, 3.9, 4.3, 4.1, 4.0];
delV_WC = [5.5, 4.6, 3.9, 3.2, 2.6, 2.0, 1.8, 1.3, 1.0, 0.9, 0.6];

% Initial Conc of NaOH
CA0 = 3.2; % mol/L

% Initial Conc of CH3COO-C4H9 (aq)
CB0_1 = 0;

% Initial Conc of CH3COO-C4H9 (organic)
CB0_2 = 7.6; % mol/L

% Initial Conc vector
C0 = [CA0, CB0_1, CB0_2];

%-------------------------------------------------------------------------%

% VOLUME -> CONCENTRATION (CA)
CA_NC = 0.5 .* delV_NC;
CA_WC = 0.5 .* delV_WC;

% CONCENTRATION (CA) -> CONVERSION (XA)
XA_NC = 1 - (CA_NC ./ CA0);
XA_WC = 1 - (CA_WC ./ CA0);

%-------------------------------------------------------------------------%

% PLOTTING

% concentration vs time plot
figure(1)
hold on
plot(t, CA_NC, ...
    Color="Red", LineWidth=1.75, Marker="o", ...
    DisplayName="No Catalyst")
plot(t, CA_WC, ...
    Color="Blue", LineWidth=1.75, Marker="o", ...
    DisplayName="With Catalyst")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Concentration (mol/L)')
fontsize(20, "points")

% conversion vs time plot
figure(2)
hold on
plot(t, XA_NC, ...
    Color="Red", LineWidth=1.75, Marker="o", ...
    DisplayName="No Catalyst")
plot(t, XA_WC, ...
    Color="Blue", LineWidth=1.75, Marker="o", ...
    DisplayName="With Catalyst")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Conversion')
fontsize(20, "points")

%-------------------------------------------------------------------------%

% CHOOSE CASE (CONCENTRATION)

% CA = CA_NC; % no catalyst
CA = CA_WC; % with catalyst

%-------------------------------------------------------------------------%

% DIFFERENTIAL ANALYSIS

% derivative of concentration in (mol/Lmin)
dCAdt = zeros(size(CA));

% forward diff
dCAdt(1) = (CA(2) - CA(1)) / (t(2) - t(1));

% central diff
for i = 2:length(t)-1
    dCAdt(i) = (CA(i+1) - CA(i-1)) / (t(i+1) - t(i-1));
    if dCAdt(i) >= 0
        dCAdt(i) = (CA(i+1) - CA(i)) / (t(i+1) - t(i));
    end
    if dCAdt(i) >= 0
        dCAdt(i) = (CA(i) - CA(i-1)) / (t(i) - t(i-1));
    end
end

% backward diff
dCAdt(end) = (CA(end) - CA(end-1)) / (t(end) - t(end-1));

% fitting polynomial
xdata = log10(CA);
ydata = log10(-dCAdt);
DA_p = polyfit(xdata, ydata, 1);

% plotting derivative vs conc in log-log
x = linspace(-0.6, 0.6, 100); % with catalyst
% x = linspace(0.2, 0.5, 100); % no catalyst
figure(3)
hold on
scatter(xdata, ydata, 50, "red", "filled", DisplayName="Data Points")
plot(x, polyval(DA_p, x), "Color", "blue", "LineWidth", 1.75, DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel('Concentration (mol/L)')
ylabel('-1 * Derivative (mol/L*min)')
fontsize(20, "points")

% Parameter values
DA_a = DA_p(1);
DA_kr = 10^(DA_p(2));

% Goodness of fit
DA_r2 = goodnessFit(ydata, polyval(DA_p, xdata));

%-------------------------------------------------------------------------%

% INTEGRAL ANALYSIS

INT_CA = zeros(size(CA));

for i = 2:length(t)
    INT_CA(2:i) = INT_CA(2:i) + ((CA(i) .^ (-DA_a)) * (CA(i) - CA(i-1)));
end

INT_CA = INT_CA(2:end);
INT_p = polyfit(t(2:end), INT_CA, 1);

% plotting integral vs time
x = linspace(0, 55, 100);
figure(4)
hold on
scatter(t(2:end), INT_CA, 50, "red", "filled", DisplayName="Data Points")
plot(x, polyval(INT_p, x), "Color", "blue", "LineWidth", 1.75, DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Integral')
fontsize(20, "points")

INT_kr = INT_p(1);

% Goodness of fit
INT_r2 = goodnessFit(INT_CA, polyval(INT_p, t(2:end)));

%-------------------------------------------------------------------------%

% REGRESSION ANALYSIS

% parameters given in paper (taken as guess here)
p0 = [0.076, 0.0093, 2.0, 0.1];

% estimating (fitting) values using regression
options = optimoptions("lsqnonlin", ...
    "Algorithm", "levenberg-marquardt", ...
    "FunctionTolerance", 1e-8, ...
    "MaxIterations", 1e4, ...
    "Display", "none");
p = lsqnonlin(@(p) errorFunc(C0, CA, t, p), ...
    p0, [0, 0, 0, 0], [], [], [], [], [], [], options);

% solving ode
tspan = 0:5:50;
[ode_t, ode_C] = ode45(@(t, C) reaxEq(C, t, p), tspan, C0);

% plotting graph
figure(5)
hold on
plot(ode_t, ode_C(:, 1), LineWidth=1.75, Marker="o", DisplayName="CA")
plot(ode_t, ode_C(:, 2), LineWidth=1.75, Marker="o", DisplayName="CB (1)")
plot(ode_t, ode_C(:, 3), LineWidth=1.75, Marker="o", DisplayName="CB (2)")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Concentration (mol/L)')
fontsize(20, "points")

% plotting graph
figure(6)
hold on
plot(ode_t, ode_C(:, 1), LineWidth=1.75, Marker="o", DisplayName="Estimated CA")
scatter(t, CA, 50, "red", "filled", DisplayName="True CA")
hold off
grid on
legend(Location="best")
xlabel('Time (min)')
ylabel('Concentration (mol/L)')
fontsize(20, "points")

% Goodness of fit
REG_r2 = goodnessFit(CA, (ode_C(:, 1))');

%-------------------------------------------------------------------------%

% LOCAL FUNCTIONS

% Model Equations (ODEs)
function [ dydt ] = reaxEq(C, t, p)

% concentration
CA = C(1);
CB_1 = C(2);
CB_2 = C(3);

% parameters
kr = p(1);
kma = p(2);
a = p(3);
b = p(4);

% provided in paper
kd = 66.0;

% output
dydt = zeros(3, 1);
dydt(1) = -1 * kr * (CA ^ a) * (CB_1 ^ b);
dydt(2) = (kr * (CA ^ a) * (CB_1 ^ b)) + ...
    (kma * ( (CB_2 / kd) - CB_1 ));
dydt(3) = -1 * kma * ( (CB_2 / kd) - CB_1 );

end

% Minimizing Function 
function [ RMSE ] = errorFunc(C0, CA, t, p)

% solving ode
[~, ode_C] = ode45(@(t, C) reaxEq(C, t, p), t, C0);
CA_hat = ode_C(:, 1);

% error
err = CA - CA_hat;

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