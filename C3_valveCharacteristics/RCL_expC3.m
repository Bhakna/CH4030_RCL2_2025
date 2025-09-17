clc;
clear variables;
close all;

% All data is of the form:
%   Column 1 = "Pressure (psi)"
%   Column 2 = "Stem Opening (units)"
%   Column 3 = "Flow Rate (LPH)"

%-------------------------------------------------------------------------%

% Data Initialization

% File Names
fileNames = ["RCL2_expC3_linearForward.xlsx"; ...
    "RCL2_expC3_linearBackward.xlsx"; ...
    "RCL2_expC3_equalForward.xlsx"; ...
    "RCL2_expC3_equalBackward.xlsx"; ...
    "RCL2_expC3_quickForward.xlsx"; ...
    "RCL2_expC3_quickBackward.xlsx"];


% Max Stem Position in Units
maxUnits = 30;

% Max Flow Rate in LPH
maxFlow = 850;

% Curve Fitting Parameters
params = zeros(6, 2);
RMSE = zeros(6, 1);
R2 = zeros(6, 1);

%-------------------------------------------------------------------------%

% Curve Fitting and Figure Plotting

for i = 1:6

    % getting data from .xlsx file
    saveFile = fileNames(i);
    data = table2array(readtable(saveFile));

    % Control Valve Data Name
    valveName1 = char(saveFile);
    valveName1 = string(valveName1(12:end-5));

    % getting "Stem Position" and "Flow Rate" as x and y
    xdata1 = data(:, 2);
    ydata1 = data(:, 3);

    % Excluding x=0 points
    includePoints = [];
    for j = 1:length(xdata1)
        if xdata1(j) >= 1
            includePoints = [includePoints, j];
        end
    end

    % Cleaned and Normalized Data
    xdata1 = xdata1(includePoints) ./ maxUnits;
    ydata1 = ydata1(includePoints) ./ maxFlow;

    % Fitting Parameters
    p0 = [-1, -1];
    params(i, :) = fminunc(@(p) rmse(p, xdata1, ydata1), p0);

    % Goodness of Fit
    RMSE(i) = rmse(params(i, :), xdata1, ydata1);
    R2(i) = r2(params(i, :), xdata1, ydata1);

end

for i = 1:2:6

    % getting data from .xlsx file
    saveFile = fileNames(i);
    data = table2array(readtable(saveFile));

    % Control Valve Data Name
    valveName1 = char(saveFile);
    valveName1 = string(valveName1(12:end-5));

    % Fitted Data
    xfit1 = linspace(0, 1, 50);
    yfit1 = genFunc(params(i, :), xfit1);

    % Original Data (normalized)
    xdata1 = data(:, 2) ./ maxUnits;
    ydata1 = data(:, 3) ./ maxFlow;

    % case 2

    % getting data from .xlsx file
    saveFile = fileNames(i + 1);
    data = table2array(readtable(saveFile));

    % Control Valve Data Name
    valveName2 = char(saveFile);
    valveName2 = string(valveName2(12:end-5));

    % Fitted Data
    xfit2 = linspace(0, 1, 50);
    yfit2 = genFunc(params(i + 1, :), xfit1);

    % Original Data (normalized)
    xdata2 = data(:, 2) ./ maxUnits;
    ydata2 = data(:, 3) ./ maxFlow;

    % Figure
    figPlotter(xdata1, ydata1, xfit1, yfit1, valveName1, ...
        xdata2, ydata2, xfit2, yfit2, valveName2)

end


%-------------------------------------------------------------------------%

% % Composite Figure
% 
% close all;
% xdata = linspace(0, 1, 100);
% 
% figure(7)
% hold on
% for i = 1:2:6
%     ydata = genFunc(params(i, :), xdata);
%     plot(xdata.*100, ydata.*100, LineWidth=2, ...
%         MarkerIndices=[1, 10:10:100], Marker="o")
% end
% hold off
% grid on
% legend('Linear', 'Equal', 'Quick', Location='best')
% axis([-0.05 1.05 -0.05 1.05].*100)
% xlabel('Stem Position (%)')
% ylabel('Flow Rate (%)')
% title('Forward')
% fontsize(24, "points")
% 
% figure(8)
% hold on
% for i = 2:2:6
%     ydata = genFunc(params(i, :), xdata);
%     plot(xdata.*100, ydata.*100, LineWidth=2, ...
%         MarkerIndices=[1, 10:10:100], Marker="o")
% end
% hold off
% grid on
% legend('Linear', 'Equal', 'Quick', Location='best')
% axis([-0.05 1.05 -0.05 1.05].*100)
% xlabel('Stem Position (%)')
% ylabel('Flow Rate (%)')
% title('Backward')
% fontsize(24, "points")

%-------------------------------------------------------------------------%

% Local Functions

% general function
function [ f_x ] = genFunc(p, x)
d = -1 * p(1);
c = 1 + p(1) - (p(1) * exp(p(2)));
f_x = (p(1) .* exp(p(2) .* x)) + (c .* x) + d;
end

% minimizing function
function [ output ] = rmse(p, x, y)
y_hat = genFunc(p, x);
err = y - y_hat;
output = sqrt(mean(err.^2, "all"));
end

% calculating goodness of fit
function [ output ] = r2(p, x, y)
y_hat = genFunc(p, x);
RSS = sum((y - y_hat).^2, "all");
TSS = sum((y - mean(y, "all")).^2, "all");
output = 1 - (RSS / TSS);
end

% figure plotter
function [] = figPlotter(xdata1, ydata1, xfit1, yfit1, valveName1, ...
        xdata2, ydata2, xfit2, yfit2, valveName2)

figure("Name", valveName1)

hold on

plot(xfit1.*100, yfit1.*100, LineWidth=2, Color="Blue", DisplayName="Fitted Curve")
plot(xfit2.*100, yfit2.*100, LineWidth=2, Color="Red", DisplayName="Fitted Curve")

scatter(xdata1.*100, ydata1.*100, 50, "Blue", "filled", DisplayName=valveName1)
scatter(xdata2.*100, ydata2.*100, 50, "Red", "filled", DisplayName=valveName2)

hold off

grid on
legend(Location="best")
axis([-0.05 1.05 -0.05 1.05].*100)
xlabel('Stem Position (%)')
ylabel('Flow Rate (%)')

fontsize(24, "points")

end

%-------------------------------------------------------------------------%