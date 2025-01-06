% Sample data vector
%data = b1; % Replace this with your actual data vector


function [f, xi]=nonparamtericDistfitter(data)
    % fitNonParametricDistribution Fits a non-parametric distribution to the data
    % and optionally saves the plot to a file.
    %
    % Inputs:
    %   data - Vector of data points
    %   saveFileName - (Optional) File name to save the plot
    %
    % Example usage:
    %   data = randn(1000, 1); % Generate some random data
    %   fitNonParametricDistribution(data, 'kde_fit.png');

% Perform kernel density estimation
[f, xi] = ksdensity(data);
%%%% you can also use the following pd1 = fitdist(b1,'kernel','kernel','normal','support','unbounded');
% Plot the data and the fitted density
figure;
histogram(data, 'Normalization', 'pdf'); % Plot the histogram of the data
hold on;
plot(xi, f, 'LineWidth', 2); % Plot the fitted density
hold off;

% Add labels and title
xlabel('Data Values');
ylabel('Density');
title('Kernel Density Estimation');
legend('Data Histogram', 'Kernel Density Estimate');

% Optionally, save the plot as an image file
saveas(gcf, 'kde_fit.png');
end