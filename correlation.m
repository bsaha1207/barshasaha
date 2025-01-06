clc;
clear;

% Loop through each case file from 1 to 7
for caseNum = 1:7
    % Load the estimated parameters from the file
    caseNum
    fileName = sprintf('Mcase%d process.xls', caseNum);
    caseData = readmatrix(fileName);
    
    % Extract the estimated parameters from columns 14 to 26
    params = caseData(:, 14:26);
    
    % Define parameter names for reference in the plot titles
    paramNames = {'a', 'w', 'bs', 'ba', 'zi', 'e', 'u', 'sigma', 'r', 'eta', 'del', 'phi'};
    
    % Compute the correlation matrix for the parameters
    R = corrcoef(params);
    
    % Display the correlation matrix for each case
    fprintf('Correlation matrix for Case %d:\n', caseNum);
    disp(R);
    
    % Find parameter pairs with abs(correlation) >= 0.6
    [rowIdx, colIdx] = find(abs(R) >= 0.6 & abs(R) < 1); % Avoid self-correlation (1 on the diagonal)
    
    % Loop through each parameter pair and create scatter plots
    for k = 1:length(rowIdx)
        % Define the parameters to plot based on correlation index
        xParam = params(:, rowIdx(k));
        yParam = params(:, colIdx(k));
        
        % Create a scatter plot for the parameter pair
        figure;
        scatter(xParam, yParam, 'filled');
        
        % Label the axes and add title
        xlabel(paramNames{rowIdx(k)}, 'FontSize', 12);
        ylabel(paramNames{colIdx(k)}, 'FontSize', 12);
        title(sprintf('Case %d: Scatter Plot of %s vs %s (|corr|=%.2f)', caseNum, ...
                      paramNames{rowIdx(k)}, paramNames{colIdx(k)}, abs(R(rowIdx(k), colIdx(k)))), 'FontSize', 14);
        grid on;
        
        % Calculate and display p-value for this correlation
        [~, pValue] = corr(xParam, yParam);
        fprintf('Case %d: p-value for %s vs %s = %.4f\n', caseNum, ...
                paramNames{rowIdx(k)}, paramNames{colIdx(k)}, pValue);
        
        % Optionally, save each scatter plot (e.g., as PNG)
        saveas(gcf, sprintf('Case%d_Scatter_%s_vs_%s.png', caseNum, paramNames{rowIdx(k)}, paramNames{colIdx(k)}));
    end
end
