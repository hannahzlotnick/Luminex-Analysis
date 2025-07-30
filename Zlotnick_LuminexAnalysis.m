%%%%%%% Part 1: Fit the Curve %%%%%%%

%Key%
% MIP1a 1-2, IL1B 3-4, IL4 5-6, IP10 7-8, IL6 9-10, IL8 11-12, IL10 13-14, 
% IL12 15-16, IL13 17-18, IL17 19-20, IFNG 21-22, GMCSF 23-24, TNFa 25-26, 
% MIP1B 27-28, IFNa 29-30, MCP1 31-32, Pselectin 33-34, IL1a 35-36, 
% ICAM1 37-38, Eselectin 39-40.

standards = 'Standard curves.xlsx';
allstd = readmatrix(standards);
x_data = allstd(:,13); %Should be the known concentrations in pg mL (independent variable)
y_data = allstd(:,14); %Median values from luminex (dependent variable)

% Initial parameter estimates
A_init = min(y_data);
D_init = max(y_data);
C_init = mean(x_data);
B_init = 1;

initialParams = [A_init, B_init, C_init, D_init];

% Fit the model
[paramFit, R, J, CovB, MSE] = nlinfit(x_data, y_data, ...
    @(beta, x) logistic4pl(x, beta(1), beta(2), beta(3), beta(4)), initialParams);

% Generate fitted values
y_fit = logistic4pl(x_data, paramFit(1), paramFit(2), paramFit(3), paramFit(4));

% Calculate residuals
residuals = y_data - y_fit;

% Calculate the Total Sum of Squares (SST)
SST = sum((y_data - mean(y_data)).^2);

% Calculate the Residual Sum of Squares (RSS)
RSS = sum(residuals.^2);

% Calculate R^2
R_squared = 1 - (RSS / SST);

% Display R^2
disp(['R^2 = ', num2str(R_squared)]);

% Plot the data and the fitted curve
figure;
scatter(x_data, y_data, 'filled'); % Original data
hold on;
plot(x_data, y_fit, 'r-', 'LineWidth', 2); % Fitted curve
set(gca, 'XScale', 'log'); % Log scale for better visualization
xlabel('X (log scale)');
ylabel('Y');
title('4-Parameter Logistic Curve Fit');
legend('Data', 'Fitted Curve');
hold off;

%%%%%%%%%%%%%%%%%% Part 2: Get the Concentration Data from Median Values%%%%%%%%%%%%%%%%%%%%

% Parameters from the fitted model
A_fit = paramFit(1);
B_fit = paramFit(2);
C_fit = paramFit(3);
D_fit = paramFit(4);

% Given list of median values
filename = 'Median Values_2.xlsx';
data = readmatrix(filename);
%Column 2: MIP1a
%Column 4: IL-1B
%Column 6: IL-4
%Column 8: IP-10
%Column 10: IL-6
%Column 12: IL-8
%Column 14: IL-10
%Column 16: IL-12p70
%Column 18: IL-13
%Column 20: IL-17a
%Column 22: IFNg
%Column 24: GMCSF
%Column 26: TNFa
%Column 28: MIP1b
%Column 30: IFNa
%Column 32: MCP1
%Column 34: P-selectin
%Column 36: IL-1a
%Column 38: ICAM-1
%Column 40: E-selectin
column = data(:,14); %Change column number accordingly
median_values = column;

% Preallocate x-values array
x_values = zeros(size(median_values));

% Calculate corresponding x-values
% Calculate corresponding x-values
for i = 1:length(median_values)
    y_f = median_values(i);
    
    % Check if the y-value is within the valid range
    if y_f < A_fit || y_f > D_fit
        warning('y-value %f is out of the valid range [%f, %f].', y_f, A_fit, D_fit);
        x_values(i) = NaN; % Assign NaN for out-of-range values
        continue; % Skip to the next iteration
    end
    
    % Calculate corresponding x
    x_values(i) = C_fit * ((y_f - A_fit) / (D_fit - y_f))^(1 / B_fit);
end

% Display the computed x-values
disp('Computed x-values:');
disp(x_values);

