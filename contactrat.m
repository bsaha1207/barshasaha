% Given data points for percentage avoiding public places (approximate values)
months = datetime(2020, 4, 1) + calmonths(0:8); % From Apr 2020 to Jan 2022
percent_avoiding = [79, 53, 61, 44, 27, 35, 29, 41, 28]; % Actual given values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 4, 1) + calmonths(0:35); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the existing data to extrapolate
p = polyfit(days_since_start, percent_avoiding, 2); % Using a 2nd-degree polynomial

% Generate extended data using the polynomial fit
extended_percent_avoiding = polyval(p, extended_days_since_start);

% Ensure the extrapolated values follow a decline
extended_percent_avoiding(extended_percent_avoiding < 0) = 0; % Prevent negative values

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extended_percent_avoiding, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values
%%
% Given data points for percentage avoiding public places (detailed values from the graph)
months = datetime(2020, 1, 1) + calmonths(0:23); % From Jan 2020 to Jan 2022
percent_avoiding = [23, 70, 79, 65, 60, 58, 56, 54, 53, 51, 50, 57, 55, 50, 45, 42, 40, 37, 35, 33, 29, 31, 34, 28]; % Detailed values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 1, 1) + calmonths(0:51); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the observed data
p = polyfit(days_since_start, percent_avoiding, 4); % Using a 4th-degree polynomial to capture the trend

% Generate data using the polynomial fit for the entire period
extended_percent_avoiding = polyval(p, extended_days_since_start);

% Manually adjust the trend to avoid unrealistic values after January 2022
for i = length(days_since_start)+1:length(extended_days_since_start)
    if extended_percent_avoiding(i) < 20
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) - 0.2; % Gradual decline
    end
    if extended_percent_avoiding(i) < 0
        extended_percent_avoiding(i) = 0; % Ensure it doesn't go below 0
    end
end

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extended_percent_avoiding, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values

%%
% Given data points for percentage avoiding public places (detailed values from the graph)
months = datetime(2020, 1, 1) + calmonths(0:23); % From Jan 2020 to Jan 2022
percent_avoiding = [23, 70, 79, 65, 60, 58, 56, 54, 53, 51, 50, 57, 55, 50, 45, 42, 40, 37, 35, 33, 29, 31, 34, 28]; % Detailed values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 1, 1) + calmonths(0:51); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the observed data
p = polyfit(days_since_start, percent_avoiding, 4); % Using a 4th-degree polynomial to capture the trend

% Generate data using the polynomial fit for the entire period
extended_percent_avoiding = polyval(p, extended_days_since_start);

% Simulate the effect of the Omicron variant and other fluctuations
for i = length(days_since_start)+1:length(extended_days_since_start)
    if extended_days_since_start(i) >= days(datetime(2022, 12, 1) - months(1)) && extended_days_since_start(i) <= days(datetime(2023, 3, 1) - months(1))
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) + 5; % Simulate increase during Omicron
    elseif extended_percent_avoiding(i) < extended_percent_avoiding(i-1) && extended_percent_avoiding(i) < 20
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) - 0.1; % Gradual decline after Omicron
    end
    if extended_percent_avoiding(i) < 0
        extended_percent_avoiding(i) = 0; % Ensure it doesn't go below 0
    end
    if extended_percent_avoiding(i) > 100
        extended_percent_avoiding(i) = 100; % Ensure it doesn't exceed 100%
    end
end

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extended_percent_avoiding, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values
%%
% Given data points for percentage avoiding public places (detailed values from the graph)
months = datetime(2020, 1, 1) + calmonths(0:23); % From Jan 2020 to Jan 2022
percent_avoiding = [23, 70, 79, 65, 60, 58, 56, 54, 53, 51, 50, 57, 55, 50, 45, 42, 40, 37, 35, 33, 29, 31, 34, 28]; % Detailed values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 1, 1) + calmonths(0:51); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the observed data
p = polyfit(days_since_start, percent_avoiding, 4); % Using a 4th-degree polynomial to capture the trend

% Generate data using the polynomial fit for the entire period
extended_percent_avoiding = polyval(p, extended_days_since_start);

% Simulate the effect of the Omicron variant and other fluctuations
for i = length(days_since_start)+1:length(extended_days_since_start)
    % Introduce a gradual increase during the first constant part
    if extended_days_since_start(i) > days(datetime(2022, 4, 1) - months(1)) && extended_days_since_start(i) <= days(datetime(2022, 12, 1) - months(1))
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) + 0.2; % Gradual increase after the dip
    elseif extended_days_since_start(i) > days(datetime(2022, 12, 1) - months(1)) && extended_days_since_start(i) <= days(datetime(2023, 3, 1) - months(1))
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) + 5; % Simulate increase during Omicron
    elseif extended_days_since_start(i) > days(datetime(2023, 3, 1) - months(1)) && extended_days_since_start(i) <= days(datetime(2023, 8, 1) - months(1))
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) + 0.3; % Gradual increase after Omicron
    elseif extended_days_since_start(i) > days(datetime(2023, 8, 1) - months(1))
        extended_percent_avoiding(i) = extended_percent_avoiding(i-1) - 0.4; % Gradual decline after the increase
    end
    
    % Ensure values stay within realistic bounds and do not go to zero prematurely
    if extended_percent_avoiding(i) < 10
        extended_percent_avoiding(i) = 10; % Prevent values from dropping below 10%
    end
    if extended_percent_avoiding(i) > 100
        extended_percent_avoiding(i) = 100; % Prevent values from exceeding 100%
    end
end

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extended_percent_avoiding, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values
%%
% Given data points for percentage avoiding public places (detailed values from the graph)
months = datetime(2020, 1, 1) + calmonths(0:23); % From Jan 2020 to Jan 2022
percent_avoiding = [23, 70, 79, 65, 60, 58, 56, 54, 53, 51, 50, 57, 55, 50, 45, 42, 40, 37, 35, 33, 29, 31, 34, 28]; % Detailed values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 1, 1) + calmonths(0:51); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the observed data
p = polyfit(days_since_start, percent_avoiding, 3); % Using a 3rd-degree polynomial to capture the trend

% Generate data using the polynomial fit for the observed period
observed_fit = polyval(p, days_since_start);

% Extrapolate data after January 2022
extrapolated_data = polyval(p, extended_days_since_start);

% Modify the extrapolation to reflect a gradual increase after January 2022
for i = length(days_since_start)+1:length(extended_days_since_start)
    % Simulate an increase after January 2022
    extrapolated_data(i) = extrapolated_data(i-1) - 0.5; % Gradual decrease to fit the trend of not avoiding
    if extrapolated_data(i) < 0
        extrapolated_data(i) = 0; % Prevent it from going negative
    end
end

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extrapolated_data, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values
%%
% Given data points for percentage avoiding public places (detailed values from the graph)
months = datetime(2020, 1, 1) + calmonths(0:23); % From Jan 2020 to Jan 2022
percent_avoiding = [23, 70, 79, 65, 60, 58, 56, 54, 53, 51, 50, 57, 55, 50, 45, 42, 40, 37, 35, 33, 29, 31, 34, 28]; % Detailed values

% Create a time vector extending to March 2024
extended_months = datetime(2020, 1, 1) + calmonths(0:51); % Extend to Mar 2024

% Calculate the number of days since the first month for polyfit
days_since_start = days(months - months(1)); % Days since the start of the data
extended_days_since_start = days(extended_months - months(1)); % Extended days

% Fit a polynomial to the observed data
p = polyfit(days_since_start, percent_avoiding, 3); % Using a 3rd-degree polynomial to capture the trend

% Generate data using the polynomial fit for the observed period
observed_fit = polyval(p, days_since_start);

% Extrapolate data after January 2022
extrapolated_data = polyval(p, extended_days_since_start);

% Modify the extrapolation to reflect a gradual increase after January 2022
for i = length(days_since_start)+1:length(extended_days_since_start)
    % Simulate an increase after January 2022
    extrapolated_data(i) = extrapolated_data(i-1) + 0.5; % Gradual increase after January 2022
    if extrapolated_data(i) > 100
        extrapolated_data(i) = 100; % Ensure it doesn't exceed 100%
    end
end

% Plot the original data
figure;
plot(months, percent_avoiding, 'g', 'LineWidth', 2); hold on;

% Plot the extended data
plot(extended_months, extrapolated_data, 'g--', 'LineWidth', 2);

% Annotate the plot
xlabel('Date');
ylabel('% Avoiding Public Places');
title('Americans Social Distancing Practices During the Pandemic');
legend('Observed Data', 'Extrapolated Data');
grid on;

% Set the x-axis to show every 3 months for clarity
xtickformat('MMM yyyy');
xticks(extended_months(1):calmonths(3):extended_months(end));

% Display the plot
xlim([months(1) extended_months(end)]);
ylim([0 100]); % Set y-axis limits to avoid unrealistic values
%%
clc;
close all;

% Given values for ca for each case
caValues = [0.21, 0.47, 0.6, 0.65, 0.58, 0.54, 0.4];

% Define indices for bs and ba in the files
bsIndex = 17;
baIndex = 18;

% Initialize arrays to store pa and ps values for all cases
paValues = {};
psValues = {};

% Loop through each case to calculate pa and ps
for caseNum = 1:numCases
    % Construct the file name for input
    inputFile = sprintf('Mcase%d process.xls', caseNum);
    
    % Read the data from the existing Excel file
    Mestimate = readmatrix(inputFile);
    
    % Extract bs and ba data
    bsData = Mestimate(:, bsIndex);
    baData = Mestimate(:, baIndex);
    
    % Calculate pa and ps
    ca = caValues(caseNum);
    pa = baData / ca;
    ps = 2 * bsData / ca;
    
    % Store the calculated pa and ps values
    paValues{caseNum} = pa;
    psValues{caseNum} = ps;
    
    % Save each case's pa and ps to a separate sheet in the same Excel file
    writematrix(pa, 'paValues.xlsx', 'Sheet', sprintf('Case %d', caseNum));
    writematrix(ps, 'psValues.xlsx', 'Sheet', sprintf('Case %d', caseNum));
end

% Combine all cases into one array for boxplot and histogram
allPaData = vertcat(paValues{:});
allPsData = vertcat(psValues{:});

% Generate box plots for pa and ps
figure;
boxplot(allPaData, 'Labels', arrayfun(@(x) sprintf('Case %d', x), 1:numCases, 'UniformOutput', false));
title('Boxplot of pa Across All Cases');
xlabel('Cases');
ylabel('pa');

figure;
boxplot(allPsData, 'Labels', arrayfun(@(x) sprintf('Case %d', x), 1:numCases, 'UniformOutput', false));
title('Boxplot of ps Across All Cases');
xlabel('Cases');
ylabel('ps');

% Generate histograms for pa and ps with uniform bin edges
figure;
binEdges = 0:0.05:1; % Adjust the bin width and range as needed
histogram(allPaData, 'Normalization', 'probability', 'BinEdges', binEdges);
title('Normalized Histogram of pa Across All Cases');
xlabel('pa');
ylabel('Probability');
grid on;

figure;
histogram(allPsData, 'Normalization', 'probability', 'BinEdges', binEdges);
title('Normalized Histogram of ps Across All Cases');
xlabel('ps');
ylabel('Probability');
grid on;
%%
clc;
close all;

% Given values for ca for each case
caValues = [0.21, 0.47, 0.6, 0.65, 0.58, 0.54, 0.4];

% Define indices for bs and ba in the files
bsIndex = 17;
baIndex = 18;

% Initialize arrays to store pa and ps values for all cases
paValues = {};
psValues = {};

% Initialize arrays to store grouping information
paGroups = [];
psGroups = [];

% Loop through each case to calculate pa and ps
for caseNum = 1:numCases
    % Construct the file name for input
    inputFile = sprintf('Mcase%d process.xls', caseNum);
    
    % Read the data from the existing Excel file
    Mestimate = readmatrix(inputFile);
    
    % Extract bs and ba data
    bsData = Mestimate(:, bsIndex);
    baData = Mestimate(:, baIndex);
    
    % Calculate pa and ps
    ca = caValues(caseNum);
    pa = baData / ca;
    ps = 2 * bsData / ca;
    
    % Store the calculated pa and ps values
    paValues{caseNum} = pa;
    psValues{caseNum} = ps;
    
    % Save each case's pa and ps to a separate sheet in the same Excel file
    writematrix(pa, 'paValues.xlsx', 'Sheet', sprintf('Case %d', caseNum));
    writematrix(ps, 'psValues.xlsx', 'Sheet', sprintf('Case %d', caseNum));
    
    % Append the case number to the grouping variable
    paGroups = [paGroups; repmat(caseNum, size(pa))];
    psGroups = [psGroups; repmat(caseNum, size(ps))];
end

% Combine all cases into one array for boxplot and histogram
allPaData = vertcat(paValues{:});
allPsData = vertcat(psValues{:});

% Generate box plots for pa and ps
figure;
boxplot(allPaData, paGroups, 'Labels', arrayfun(@(x) sprintf('Case %d', x), 1:numCases, 'UniformOutput', false));
title('Boxplot of pa Across All Cases');
xlabel('Cases');
ylabel('pa');

figure;
boxplot(allPsData, psGroups, 'Labels', arrayfun(@(x) sprintf('Case %d', x), 1:numCases, 'UniformOutput', false));
title('Boxplot of ps Across All Cases');
xlabel('Cases');
ylabel('ps');

% Generate histograms for pa and ps with uniform bin edges
figure;
binEdges = 0:0.05:1; % Adjust the bin width and range as needed
histogram(allPaData, 'Normalization', 'probability', 'BinEdges', binEdges);
title('Normalized Histogram of pa Across All Cases');
xlabel('pa');
ylabel('Probability');
grid on;

figure;
histogram(allPsData, 'Normalization', 'probability', 'BinEdges', binEdges);
title('Normalized Histogram of ps Across All Cases');
xlabel('ps');
ylabel('Probability');
grid on;
%%
clc;
%close all;

% Given values for ca for each case
csValues = [4, 2, 4.35, 5.3, 5.45, 5.25, 5.6];

% Define indices for bs and ba in the files
bsIndex = 17;

% Loop through each case to calculate ps and save them in separate Excel files
for caseNum = 1:7
    % Construct the file name for input
    inputFile = sprintf('Mcase%d process.xls', caseNum);
    
    % Read the data from the existing Excel file
    Mestimate = readmatrix(inputFile);
    
    % Extract bs data
    bsData = Mestimate(:, bsIndex);
    
    % Calculate ps
    cs = csValues(caseNum);
    ps =  bsData / cs;
    
    % Save the ps values to a separate Excel file for each case
    filename = sprintf('ps_values_%d.xlsx', caseNum);
    writematrix(ps, filename);
end

% Initialize a cell array to hold the ps values for each case
ps_values = cell(7, 1);

% Loop through each case to load the ps data
for caseNum = 1:7
    % Load the ps values from the corresponding Excel file
    filename = sprintf('ps_values_%d.xlsx', caseNum);
    ps_values{caseNum} = readmatrix(filename);
end

% Create a figure for the boxplots
figure;
hold on;

% Plot each case's boxplot with an offset
positions = 1:7; % Positions for each boxplot
for caseNum = 1:7
    % Plot the boxplot for each case
    boxplot(ps_values{caseNum}, 'Positions', positions(caseNum), 'Colors', 'b');
end

% Customize the x-axis labels and add title and labels
set(gca, 'XTick', 1:7, 'XTickLabel', {'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5', 'Case 6', 'Case 7'});
ylim([0 0.2])
xlabel('Case Number');
ylabel('ps Value');
title('Boxplot of ps Values for Different COVID-19 Variants');
grid on;
hold off;


% Define custom bin edges between 0 and 1
binEdges = 0:0.012:0.3; % Adjust the bin width (0.012) as needed

figure; % Create a new figure window

% Loop through each case
for caseNum = 1:7
    % Load the ps values from the corresponding Excel file
    filename = sprintf('ps_values_%d.xlsx', caseNum);
    psData = readmatrix(filename);
    
    % Plot the histogram with normalization and custom bin edges
    h = histogram(psData, 'Normalization', 'probability', 'BinEdges', binEdges);
    hold on; % Hold on to plot multiple histograms in the same figure
end

% Add title and labels
title('Normalized Histogram of ps Across All Cases', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('ps', 'FontSize', 11);
ylabel('Probability', 'FontSize', 11);

% Add grid for better visualization
grid on;

% Add a legend to identify each case
legend(arrayfun(@(x) sprintf('Case %d', x), 1:7, 'UniformOutput', false), 'Location', 'best');

%%
clc;
%close all;

% Given values for ca for each case
caValues = [7.5, 8.34, 7.87, 9.86, 10.09, 9.65, 10.12];
%caValues=[6.825 15.257 12.0715 16.9 18.85 23.4 26]

% Define indices for bs and ba in the files
baIndex = 18;

% Loop through each case to calculate ps and save them in separate Excel files
for caseNum = 1:7
    % Construct the file name for input
    inputFile = sprintf('Mcase%d process.xls', caseNum);
    
    % Read the data from the existing Excel file
    Mestimate = readmatrix(inputFile);
    
    % Extract bs data
    baData = Mestimate(:, baIndex);
    
    % Calculate ps
    ca = caValues(caseNum);
    pa = baData / ca;
    
    % Save the ps values to a separate Excel file for each case
    filename = sprintf('pa_values_%d.xlsx', caseNum);
    writematrix(pa, filename);
end

% Initialize a cell array to hold the ps values for each case
pa_values = cell(7, 1);

% Loop through each case to load the ps data
for caseNum = 1:7
    % Load the ps values from the corresponding Excel file
    filename = sprintf('pa_values_%d.xlsx', caseNum);
    pa_values{caseNum} = readmatrix(filename);
end

% Create a figure for the boxplots
figure;
hold on;

% Plot each case's boxplot with an offset
positions = 1:7; % Positions for each boxplot
for caseNum = 1:7
    % Plot the boxplot for each case
    boxplot(pa_values{caseNum}, 'Positions', positions(caseNum), 'Colors', 'b');
end

% Customize the x-axis labels and add title and labels
set(gca, 'XTick', 1:7, 'XTickLabel', {'Case 1', 'Case 2', 'Case 3', 'Case 4', 'Case 5', 'Case 6', 'Case 7'});
ylim([0 0.1])
xlabel('Case Number');
ylabel('pa Value');
title('Boxplot of pa Values for Different COVID-19 Variants');
grid on;
hold off;


% Define custom bin edges between 0 and 1
binEdges = 0:0.012:0.2; % Adjust the bin width (0.012) as needed

figure; % Create a new figure window

% Loop through each case
for caseNum = 1:7
    % Load the ps values from the corresponding Excel file
    filename = sprintf('pa_values_%d.xlsx', caseNum);
    paData = readmatrix(filename);
    
    % Plot the histogram with normalization and custom bin edges
    h = histogram(paData, 'Normalization', 'probability', 'BinEdges', binEdges);
    hold on; % Hold on to plot multiple histograms in the same figure
end

% Add title and labels
title('Normalized Histogram of pa Across All Cases', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('pa', 'FontSize', 11);
ylabel('Probability', 'FontSize', 11);

% Add grid for better visualization
grid on;

% Add a legend to identify each case
legend(arrayfun(@(x) sprintf('Case %d', x), 1:7, 'UniformOutput', false), 'Location', 'best');

