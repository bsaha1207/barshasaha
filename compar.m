clc;
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'bs' parameter in the file
bsIndex = 17;  % Assuming 'bs' is in the 17th column

% Define the number of cases
numCases = 7;

% Loop through each successive pair of cases (Case 1 to 2, Case 2 to 3, ..., Case 6 to 7)
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'bs' parameter data from the current and next cases
    bsData1 = data1(:, bsIndex);
    bsData2 = data2(:, bsIndex);
    
    % Remove NaN or Inf values to ensure clean data
    bsData1 = bsData1(~isnan(bsData1) & ~isinf(bsData1));
    bsData2 = bsData2(~isnan(bsData2) & ~isinf(bsData2));

    % Check the percentage of identical values
    minLength = min(length(bsData1), length(bsData2));
    numIdentical = sum(bsData1(1:minLength) == bsData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(bsData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(bsData2))]);

    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(bsData1, bsData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(bsData1, bsData2);
    
    % Display the result of the U Test
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(bsData1) > median(bsData2)
            disp('The current case has a significantly larger bs value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger bs value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in bs values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'ba' parameter in the file
baIndex = 18;  % Assuming 'ba' is in the 18th column

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'ba' parameter data
    baData1 = data1(:, baIndex);
    baData2 = data2(:, baIndex);
    
    % Remove NaN or Inf values
    baData1 = baData1(~isnan(baData1) & ~isinf(baData1));
    baData2 = baData2(~isnan(baData2) & ~isinf(baData2));

    % Check the percentage of identical values
    minLength = min(length(baData1), length(baData2));
    numIdentical = sum(baData1(1:minLength) == baData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(baData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(baData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(baData1, baData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(baData1, baData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(baData1) > median(baData2)
            disp('The current case has a significantly larger ba value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger ba value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in ba values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'zi' parameter in the file
ziIndex = 19;  % Assuming 'zi' is in the 19th column

% Define the number of cases
numCases = 7;

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'zi' parameter data
    ziData1 = data1(:, ziIndex);
    ziData2 = data2(:, ziIndex);
    
    % Remove NaN or Inf values
    ziData1 = ziData1(~isnan(ziData1) & ~isinf(ziData1));
    ziData2 = ziData2(~isnan(ziData2) & ~isinf(ziData2));

    % Check the percentage of identical values
    minLength = min(length(ziData1), length(ziData2));
    numIdentical = sum(ziData1(1:minLength) == ziData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(ziData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(ziData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(ziData1, ziData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(ziData1, ziData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(ziData1) > median(ziData2)
            disp('The current case has a significantly larger zi value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger zi value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in zi values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'e' parameter in the file
eIndex = 20;  % Assuming 'e' is in the 20th column

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'e' parameter data
    eData1 = data1(:, eIndex);
    eData2 = data2(:, eIndex);
    
    % Remove NaN or Inf values
    eData1 = eData1(~isnan(eData1) & ~isinf(eData1));
    eData2 = eData2(~isnan(eData2) & ~isinf(eData2));

    % Check the percentage of identical values
    minLength = min(length(eData1), length(eData2));
    numIdentical = sum(eData1(1:minLength) == eData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(eData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(eData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(eData1, eData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(eData1, eData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(eData1) > median(eData2)
            disp('The current case has a significantly larger e value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger e value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in e values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'eta' parameter in the file
etaIndex = 24;  % Assuming 'eta' is in the 24th column

% Define the number of cases
numCases = 7;

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'eta' parameter data
    etaData1 = data1(:, etaIndex);
    etaData2 = data2(:, etaIndex);
    
    % Remove NaN or Inf values
    etaData1 = etaData1(~isnan(etaData1) & ~isinf(etaData1));
    etaData2 = etaData2(~isnan(etaData2) & ~isinf(etaData2));

    % Check the percentage of identical values
    minLength = min(length(etaData1), length(etaData2));
    numIdentical = sum(etaData1(1:minLength) == etaData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(etaData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(etaData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(etaData1, etaData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(etaData1, etaData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(etaData1) > median(etaData2)
            disp('The current case has a significantly larger eta value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger eta value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in eta values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'delta' parameter in the file
deltaIndex = 25;  % Assuming 'delta' is in the 25th column

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'delta' parameter data
    deltaData1 = data1(:, deltaIndex);
    deltaData2 = data2(:, deltaIndex);
    
    % Remove NaN or Inf values
    deltaData1 = deltaData1(~isnan(deltaData1) & ~isinf(deltaData1));
    deltaData2 = deltaData2(~isnan(deltaData2) & ~isinf(deltaData2));

    % Check the percentage of identical values
    minLength = min(length(deltaData1), length(deltaData2));
    numIdentical = sum(deltaData1(1:minLength) == deltaData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(deltaData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(deltaData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(deltaData1, deltaData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(deltaData1, deltaData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(deltaData1) > median(deltaData2)
            disp('The current case has a significantly larger delta value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger delta value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in delta values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'phi' parameter in the file
phiIndex = 26;  % Assuming 'phi' is in the 26th column

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'phi' parameter data
    phiData1 = data1(:, phiIndex);
    phiData2 = data2(:, phiIndex);
    
    % Remove NaN or Inf values
    phiData1 = phiData1(~isnan(phiData1) & ~isinf(phiData1));
    phiData2 = phiData2(~isnan(phiData2) & ~isinf(phiData2));

    % Check the percentage of identical values
    minLength = min(length(phiData1), length(phiData2));
    numIdentical = sum(phiData1(1:minLength) == phiData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(phiData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(phiData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(phiData1, phiData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(phiData1, phiData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(phiData1) > median(phiData2)
            disp('The current case has a significantly larger phi value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger phi value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in phi values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    psData1 = readmatrix(sprintf('ps_values_%d.xlsx', caseNum));     % Read ps data for current case
    psData2 = readmatrix(sprintf('ps_values_%d.xlsx', caseNum + 1)); % Read ps data for next case
    
    % Remove NaN or Inf values to ensure clean data
    psData1 = psData1(~isnan(psData1) & ~isinf(psData1));
    psData2 = psData2(~isnan(psData2) & ~isinf(psData2));

    % Check the percentage of identical values
    minLength = min(length(psData1), length(psData2));
    numIdentical = sum(psData1(1:minLength) == psData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(psData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(psData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(psData1, psData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(psData1, psData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(psData1) > median(psData2)
            disp('The current case has a significantly larger ps value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger ps value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in ps values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    paData1 = readmatrix(sprintf('pa_values_%d.xlsx', caseNum));     % Read pa data for current case
    paData2 = readmatrix(sprintf('pa_values_%d.xlsx', caseNum + 1)); % Read pa data for next case
    
    % Remove NaN or Inf values to ensure clean data
    paData1 = paData1(~isnan(paData1) & ~isinf(paData1));
    paData2 = paData2(~isnan(paData2) & ~isinf(paData2));

    % Check the percentage of identical values
    minLength = min(length(paData1), length(paData2));
    numIdentical = sum(paData1(1:minLength) == paData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(paData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(paData2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(paData1, paData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(paData1, paData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(paData1) > median(paData2)
            disp('The current case has a significantly larger pa value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger pa value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in pa values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'R0' parameter in the file
R0Index = 27;  % Assuming 'R0' is in the 27th column

% Define the number of cases
numCases = 7;

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'R0' parameter data
    R0Data1 = data1(:, R0Index);
    R0Data2 = data2(:, R0Index);
    
    % Remove NaN or Inf values
    R0Data1 = R0Data1(~isnan(R0Data1) & ~isinf(R0Data1));
    R0Data2 = R0Data2(~isnan(R0Data2) & ~isinf(R0Data2));

    % Check the percentage of identical values
    minLength = min(length(R0Data1), length(R0Data2));
    numIdentical = sum(R0Data1(1:minLength) == R0Data2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(R0Data1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(R0Data2))]);
    
    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(R0Data1, R0Data2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(R0Data1, R0Data2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(R0Data1) > median(R0Data2)
            disp('The current case has a significantly larger R0 value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger R0 value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in R0 values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end

%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'lambda' parameter in the file
lambdaIndex = 14;  % Assuming 'lambda' is in the 14th column

% Define the number of cases
numCases = 7;


% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'lambda' parameter data
    lambdaData1 = data1(:, lambdaIndex);
    lambdaData2 = data2(:, lambdaIndex);
    
    % Remove NaN or Inf values
    lambdaData1 = lambdaData1(~isnan(lambdaData1) & ~isinf(lambdaData1));
    lambdaData2 = lambdaData2(~isnan(lambdaData2) & ~isinf(lambdaData2));

    % Check the percentage of identical values
    minLength = min(length(lambdaData1), length(lambdaData2));
    numIdentical = sum(lambdaData1(1:minLength) == lambdaData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(lambdaData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(lambdaData2))]);
    

    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(lambdaData1, lambdaData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(lambdaData1, lambdaData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(lambdaData1) > median(lambdaData2)
            disp('The current case has a significantly larger lambda value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger lambda value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in lambda values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end

%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'alpha' parameter in the file
alphaIndex = 16;  % Assuming 'alpha' is in the 16th column

% Define the number of cases
numCases = 7;

% Function to compute Cliff's Delta
% Same cliffsDelta function as before

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'alpha' parameter data
    alphaData1 = data1(:, alphaIndex);
    alphaData2 = data2(:, alphaIndex);
    
    % Remove NaN or Inf values
    alphaData1 = alphaData1(~isnan(alphaData1) & ~isinf(alphaData1));
    alphaData2 = alphaData2(~isnan(alphaData2) & ~isinf(alphaData2));

    % Check the percentage of identical values
    minLength = min(length(alphaData1), length(alphaData2));
    numIdentical = sum(alphaData1(1:minLength) == alphaData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(alphaData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(alphaData2))]);
    


    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(alphaData1, alphaData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(alphaData1, alphaData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(alphaData1) > median(alphaData2)
            disp('The current case has a significantly larger alpha value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger alpha value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in alpha values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end


%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'mu' parameter in the file
muIndex = 21;  % Assuming 'mu' is in the 21st column

% Define the number of cases
numCases = 7;

% Function to compute Cliff's Delta
% Same cliffsDelta function as before

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'mu' parameter data
    muData1 = data1(:, muIndex);
    muData2 = data2(:, muIndex);
    
    % Remove NaN or Inf values
    muData1 = muData1(~isnan(muData1) & ~isinf(muData1));
    muData2 = muData2(~isnan(muData2) & ~isinf(muData2));

    % Check the percentage of identical values
    minLength = min(length(muData1), length(muData2));
    numIdentical = sum(muData1(1:minLength) == muData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(muData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(muData2))]);
    

    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(muData1, muData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(muData1, muData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(muData1) > median(muData2)
            disp('The current case has a significantly larger mu value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger mu value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in mu values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end

%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'sigma' parameter in the file
sigmaIndex = 22;  % Assuming 'sigma' is in the 22nd column

% Define the number of cases
numCases = 7;

% Function to compute Cliff's Delta
% Same cliffsDelta function as before

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'sigma' parameter data
    sigmaData1 = data1(:, sigmaIndex);
    sigmaData2 = data2(:, sigmaIndex);
    
    % Remove NaN or Inf values
    sigmaData1 = sigmaData1(~isnan(sigmaData1) & ~isinf(sigmaData1));
    sigmaData2 = sigmaData2(~isnan(sigmaData2) & ~isinf(sigmaData2));

    % Check the percentage of identical values
    minLength = min(length(sigmaData1), length(sigmaData2));
    numIdentical = sum(sigmaData1(1:minLength) == sigmaData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(sigmaData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(sigmaData2))]);
    


    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(sigmaData1, sigmaData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(sigmaData1, sigmaData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(sigmaData1) > median(sigmaData2)
            disp('The current case has a significantly larger sigma value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger sigma value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in sigma values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end

%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'r' parameter in the file
rIndex = 23;  % Assuming 'r' is in the 23rd column

% Define the number of cases
numCases = 7;

% Function to compute Cliff's Delta
% Same cliffsDelta function as before

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'r' parameter data
    rData1 = data1(:, rIndex);
    rData2 = data2(:, rIndex);
    
    % Remove NaN or Inf values
    rData1 = rData1(~isnan(rData1) & ~isinf(rData1));
    rData2 = rData2(~isnan(rData2) & ~isinf(rData2));

    % Check the percentage of identical values
    minLength = min(length(rData1), length(rData2));
    numIdentical = sum(rData1(1:minLength) == rData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(rData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(rData2))]);
    


    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(rData1, rData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(rData1, rData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(rData1) > median(rData2)
            disp('The current case has a significantly larger r value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger r value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in r values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end
%%
close all;
clear;

% Use format longG to see full precision of p-values
format longG;

% Define the index for the 'r' parameter in the file
rIndex = 15;  % Assuming 'r' is in the 23rd column

% Define the number of cases
numCases = 7;

% Function to compute Cliff's Delta
% Same cliffsDelta function as before

% Loop through each successive pair of cases
for caseNum = 1:(numCases - 1)
    % Load the data for the current and next cases
    data1 = readmatrix(sprintf('Mcase%d process.xls', caseNum));     % Read data for current case
    data2 = readmatrix(sprintf('Mcase%d process.xls', caseNum + 1)); % Read data for next case
    
    % Extract the 'r' parameter data
    wData1 = data1(:, rIndex);
    wData2 = data2(:, rIndex);
    
    % Remove NaN or Inf values
    wData1 =wData1(~isnan(wData1) & ~isinf(wData1));
    wData2 = wData2(~isnan(wData2) & ~isinf(wData2));

    % Check the percentage of identical values
    minLength = min(length(wData1), length(wData2));
    numIdentical = sum(wData1(1:minLength) == wData2(1:minLength));
    percentageIdentical = (numIdentical / minLength) * 100;
    
    %fprintf('Percentage of identical values between Case %d and Case %d: %.2f%%\n', caseNum, caseNum + 1, percentageIdentical);

    % Display summaries for debugging
    disp(['Case ', num2str(caseNum), ' Summary:']);
    disp(['Mean: ', num2str(mean(wData1))]);
    disp(['Case ', num2str(caseNum + 1), ' Summary:']);
    disp(['Mean: ', num2str(mean(wData2))]);
    


    % Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
    [p, h, stats] = ranksum(wData1, wData2);  % ranksum for comparing two independent samples
    
    % Compute Cliff's Delta
    delta = cliffsDelta(wData1, wData2);
    
    % Display the result
    fprintf('Case %d to Case %d:\n', caseNum, caseNum + 1);
    if h == 1
        if median(wData1) > median(wData2)
            disp('The current case has a significantly larger w value than the next case (Mann-Whitney U Test).');
        else
            disp('The next case has a significantly larger w value than the current case (Mann-Whitney U Test).');
        end
    else
        disp('No significant difference in w values between the two cases (Mann-Whitney U Test).');
    end
    
    % Display the p-value
    disp(['p-value: ', num2str(p)]);
    
    % Display the rank-sum statistic
    disp(['Rank Sum Statistic: ', num2str(stats.ranksum)]);
    
    % Display Cliff's Delta
    disp(['Cliff''s Delta: ', num2str(delta)]);
    disp('-----------------------------------------------');
end