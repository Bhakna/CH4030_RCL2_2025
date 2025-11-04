clc;
close all;
clear variables;

data = readmatrix("expR4_recordedData.xlsx");
t = data(:, 1);
T = data(:, 2);

A = 39.6 * 1e-4;
T_air = 30.3;

% xdata = A .* (T - T_air);
xdata = T;
ydata = zeros(size(T));

for i = 1:length(ydata)-2
    ydata(i) = ((-1 * T(i+2)) + (4 * T(i+1)) + (-3 * T(i))) / (t(i+2) - t(i));
end
ydata(end-1) = (T(end) - T(end-1)) / (t(end) - t(end-1));
ydata(end) = ydata(end-1);

% m = sum(xdata.*ydata, "all") / sum(xdata.^2, "all");
% h = -1 * m;

fitted_param = polyfit(xdata, ydata, 1);
T_air_calc = fitted_param(2) / abs(fitted_param(1));

figure(1)
hold on
scatter(xdata, ydata, 50, "red", "filled", DisplayName="Data Points")
plot(xdata, polyval(fitted_param, xdata), ...
    "LineWidth", 1.5, "Color", "blue", DisplayName="Fitted Line")
hold off
grid on
legend(Location="best")
xlabel("Temperature")
ylabel("\partial T / \partial t")
fontsize(20, "points")

T_hat = zeros(size(T));
T_hat(1) = T(1);
for i = 2:length(T_hat)
    sTt = (fitted_param(1) * T_hat(i-1)) + fitted_param(2);
    T_hat(i) = T_hat(i-1) + (sTt * (t(i) - t(i-1)));
end

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

rhoCv = 1000 * 4180; % J/m^3 C
V = 100 * 1e-6; % m^3
H = V / A;
D = 0.155 * 1e-6;
alpha = 1/D;
limit = (1/4) * (H^2) * alpha / 10;
limit = 1 + sqrt(limit);

N = 10;
del_t = 10;
del_z = H / (N - 1);
coeff = D * del_t / (del_z ^ 2);

UA_water = abs(fitted_param(1)) * V * rhoCv; % J/s

RMSE = sqrt(mean((T - T_hat).^2, "all"));
T_mean = mean(T, "all");
R2 = 1 - (( sum((T - T_hat).^2, "all") ) / ( sum((T - T_mean).^2, "all") ));
