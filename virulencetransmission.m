% Original (t, y) data

%please udtate this code and send it back to me


%t = [1, 2, 3, 4, 5, 6, 7];  %this is the epidemic wave number
t = [92, 187, 365, 584, 765, 926, 1075]; 
%y = [2, 4, 9, 16]; % this is the mean beta_s values for each wave
load('bs_opt_values.mat')
load('ba_opt_values.mat')
load('del_opt_values.mat')
load('zi_opt_values.mat')
load('Rs_opt_values.mat')
load('R0_opt_values.mat')
load('COVIDSIRData.mat')
% Estimates
b1=bs_opt_values(:,1);
b2=bs_opt_values(:,2);
b3=bs_opt_values(:,3);
b4=bs_opt_values(:,4);
b5=bs_opt_values(:,5);
b6=bs_opt_values(:,6);
b7=bs_opt_values(:,7);


Beta_estimate1 = b1';% this is all simulates monte calro beta_s values for wave # 1
Beta_estimate2 = b2';
Beta_estimate3 = b3';
Beta_estimate4 = b4';
Beta_estimate5 = b5';
Beta_estimate6 = b6';
Beta_estimate7 = b7';


yb = [mean(Beta_estimate1), mean(Beta_estimate2),mean(Beta_estimate3),mean(Beta_estimate4),mean(Beta_estimate5),mean(Beta_estimate6),mean(Beta_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
estimates = [Beta_estimate1; Beta_estimate2; Beta_estimate3; Beta_estimate4;Beta_estimate5; Beta_estimate6;Beta_estimate7];
% Calculate the minimum and maximum bounds for shading
yb_min = min(estimates');
yb_max = max(estimates');

% Plot the original data
figure;


% Plot the shaded region
fill([t fliplr(t)], [yb_min fliplr(yb_max)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range');

% Plot the estimates
% plot(t, estimate1, '--', 'DisplayName', 'Estimate 1');
% plot(t, estimate2, '--', 'DisplayName', 'Estimate 2');
% plot(t, estimate3, '--', 'DisplayName', 'Estimate 3');
% plot(t, estimate4, '--', 'DisplayName', 'Estimate 4');

% Customize the plot
legend('show');
xlabel('days, t');
ylabel('\beta_s');
title('Original Data with Estimate Shades');
grid on;
hold on;
plot(t, yb, 'r-o', 'LineWidth', 2, 'DisplayName', 'mean \beta_s');



%%%%%%% estimates for delta

% Estimates
d1=del_opt_values(:,1);
d2=del_opt_values(:,2);
d3=del_opt_values(:,3);
d4=del_opt_values(:,4);
d5=del_opt_values(:,5);
d6=del_opt_values(:,6);
d7=del_opt_values(:,7);


Delta_estimate1 = d1';% this is all simulates monte calro delta values for wave # 1
Delta_estimate2 = d2';
Delta_estimate3 = d3';
Delta_estimate4 = d4';
Delta_estimate5 = d5';
Delta_estimate6 = d6';
Delta_estimate7 = d7';



yd = [mean(Delta_estimate1), mean(Delta_estimate2),mean(Delta_estimate3),mean(Delta_estimate4),mean(Delta_estimate5),mean(Delta_estimate6),mean(Delta_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
Destimates = [Delta_estimate1;Delta_estimate2;Delta_estimate3;Delta_estimate4;Delta_estimate5;Delta_estimate6;Delta_estimate7];

% Calculate the minimum and maximum bounds for shading
yd_min = min(Destimates');
yd_max = max(Destimates');

% Plot the original data
figure;


% Plot the shaded region
fill([t fliplr(t)], [yd_min fliplr(yd_max)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range');

% Plot the estimates
% plot(t, estimate1, '--', 'DisplayName', 'Estimate 1');
% plot(t, estimate2, '--', 'DisplayName', 'Estimate 2');
% plot(t, estimate3, '--', 'DisplayName', 'Estimate 3');
% plot(t, estimate4, '--', 'DisplayName', 'Estimate 4');

% Customize the plot
legend('show');
xlabel('days,t');
ylabel('\delta');
title('Original Data with Estimate Shades');
grid on;
hold on;
plot(t, yd, 'k-o', 'LineWidth', 2, 'DisplayName', 'mean \delta');

%%%%%%% estimates for zi

% Estimates
z1=zi_opt_values(:,1);
z2=zi_opt_values(:,2);
z3=zi_opt_values(:,3);
z4=zi_opt_values(:,4);
z5=zi_opt_values(:,5);
z6=zi_opt_values(:,6);
z7=zi_opt_values(:,7);


zi_estimate1 = z1';% this is all simulates monte calro delta values for wave # 1
zi_estimate2 = z2';
zi_estimate3 = z3';
zi_estimate4 = z4';
zi_estimate5 = z5';
zi_estimate6 = z6';
zi_estimate7 = z7';



yz = [mean(zi_estimate1), mean(zi_estimate2),mean(zi_estimate3),mean(zi_estimate4),mean(zi_estimate5),mean(zi_estimate6),mean(zi_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
zestimates = [zi_estimate1;zi_estimate2;zi_estimate3;zi_estimate4;zi_estimate5;zi_estimate6;zi_estimate7];

% Calculate the minimum and maximum bounds for shading
yz_min = min(zestimates');
yz_max = max(zestimates');

% Plot the original data
figure;


% Plot the shaded region
fill([t fliplr(t)], [yz_min fliplr(yz_max)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range');

% Plot the estimates
% plot(t, estimate1, '--', 'DisplayName', 'Estimate 1');
% plot(t, estimate2, '--', 'DisplayName', 'Estimate 2');
% plot(t, estimate3, '--', 'DisplayName', 'Estimate 3');
% plot(t, estimate4, '--', 'DisplayName', 'Estimate 4');

% Customize the plot
legend('show');
xlabel('days,t');
ylabel('\xi_{v}');
title('Original Data with Estimate Shades');
grid on;
hold on;
plot(t, yz, 'k-o', 'LineWidth', 2, 'DisplayName', 'mean \delta');


%%%%%%% estimates for Rs

% Estimates
R1=Rs_opt_values(:,1);
R2=Rs_opt_values(:,2);
R3=Rs_opt_values(:,3);
R4=Rs_opt_values(:,4);
R5=Rs_opt_values(:,5);
R6=Rs_opt_values(:,6);
R7=Rs_opt_values(:,7);


Rs_estimate1 = R1';% this is all simulates monte calro delta values for wave # 1
Rs_estimate2 = R2';
Rs_estimate3 = R3';
Rs_estimate4 = R4';
Rs_estimate5 = R5';
Rs_estimate6 = R6';
Rs_estimate7 = R7';



yr = [mean(Rs_estimate1), mean(Rs_estimate2),mean(Rs_estimate3),mean(Rs_estimate4),mean(Rs_estimate5),mean(Rs_estimate6),mean(Rs_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
Restimates = [Rs_estimate1;Rs_estimate2;Rs_estimate3;Rs_estimate4;Rs_estimate5;Rs_estimate6;Rs_estimate7];

% Calculate the minimum and maximum bounds for shading
yr_min = min(Restimates');
yr_max = max(Restimates');

% Plot the original data
figure;


% Plot the shaded region
fill([t fliplr(t)], [yr_min fliplr(yr_max)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range');

% Plot the estimates
% plot(t, estimate1, '--', 'DisplayName', 'Estimate 1');
% plot(t, estimate2, '--', 'DisplayName', 'Estimate 2');
% plot(t, estimate3, '--', 'DisplayName', 'Estimate 3');
% plot(t, estimate4, '--', 'DisplayName', 'Estimate 4');

% Customize the plot
legend('show');
xlabel('days,t');
ylabel('R0');
title('Original Data with Estimate Shades');
grid on;
hold on;
plot(t, yr, 'k-o', 'LineWidth', 2, 'DisplayName', 'mean R0');
figure
plot(t, yr, 'r-o', 'LineWidth', 2, 'DisplayName', 'mean R0');
grid on;
hold on;
plot(t, yd, 'k--o', 'LineWidth', 2, 'DisplayName', 'mean \delta');
legend('reproduction number R0s',' virulence \delta')
xlabel('days, t');
ylabel('mean values');

%%%%%%% estimates for R0

% Estimates
R01=R0_opt_values(:,1);
R02=R0_opt_values(:,2);
R03=R0_opt_values(:,3);
R04=R0_opt_values(:,4);
R05=R0_opt_values(:,5);
R06=R0_opt_values(:,6);
R07=R0_opt_values(:,7);


R0_estimate1 = R1';% this is all simulates monte calro delta values for wave # 1
R0_estimate2 = R2';
R0_estimate3 = R3';
R0_estimate4 = R4';
R0_estimate5 = R5';
R0_estimate6 = R6';
R0_estimate7 = R7';



yr0 = [mean(R0_estimate1), mean(R0_estimate2),mean(R0_estimate3),mean(R0_estimate4),mean(R0_estimate5),mean(R0_estimate6),mean(R0_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
R0estimates = [R0_estimate1;R0_estimate2;R0_estimate3;R0_estimate4;R0_estimate5;R0_estimate6;R0_estimate7];

% Calculate the minimum and maximum bounds for shading
yr0_min = min(R0estimates');
yr0_max = max(R0estimates');

% Plot the original data
figure;


% Plot the shaded region
fill([t fliplr(t)], [yr0_min fliplr(yr0_max)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range');

% Plot the estimates
% plot(t, estimate1, '--', 'DisplayName', 'Estimate 1');
% plot(t, estimate2, '--', 'DisplayName', 'Estimate 2');
% plot(t, estimate3, '--', 'DisplayName', 'Estimate 3');
% plot(t, estimate4, '--', 'DisplayName', 'Estimate 4');

% Customize the plot
legend('show');
xlabel('days,t');
ylabel('R0');
title('Original Data with Estimate Shades');
grid on;
hold on;
plot(t, yr0, 'k-o', 'LineWidth', 2, 'DisplayName', 'mean R0');
figure
plot(t, yr0, 'r-o', 'LineWidth', 2, 'DisplayName', 'mean R0');
grid on;
hold on;
plot(t, yd, 'k--o', 'LineWidth', 2, 'DisplayName', 'mean \delta');
legend('reproduction number R0',' virulence \delta')
xlabel('days, t');
ylabel('mean values');

figure
plot(t, yr0, 'r-o', 'LineWidth', 2, 'DisplayName', 'mean R0');
grid on;
hold on;
plot(t, yz, 'g--o', 'LineWidth', 2, 'DisplayName', 'mean \delta');
legend('reproduction number R0',' vacciantion \xi_{v}')
xlabel('days, t');
ylabel('mean values');

figure
sbeta=diff(yb);
sdelta=diff(yd);
diffbd=[sbeta' sdelta'];
bar(diffbd)
set(gca, 'XTickLabel', {'w1-w2', 'w2-w3', 'w3-w4', 'w4-w5', 'w5-w6','w6-w7'});
ylabel('Evolution in mean \beta and mean \delta values');
title('Relation between transmission and virulence');
legend('Change in \beta_s', 'Change in \delta', 'Location', 'Best');

figure
szi=diff(yz);
sdelta=diff(yd);
diffzd=[szi' sdelta'];
bar(diffzd)
set(gca, 'XTickLabel', {'w1-w2', 'w2-w3', 'w3-w4', 'w4-w5', 'w5-w6','w6-w7'});
ylabel('Evolution in mean \xi_{v} and mean \delta values');
title('Relation between transmission and virulence');
legend('Change in \xi_{v}', 'Change in \delta', 'Location', 'Best');


% Polynomial fitting of degree 5
p_yr = polyfit(t, yr, 5);
p_yd = polyfit(t, yd + 1.2, 6);

% Evaluate the polynomial fits
t_fit = linspace(min(t), max(t), 100);
yr_fit = polyval(p_yr, t_fit);
yd_fit = polyval(p_yd, t_fit);

% Plot the polynomial fits
figure;
plot(t, yr, 'ro', 'LineWidth', 2, 'DisplayName', 'mean R0');
hold on;
plot(t_fit, yr_fit, 'r-', 'LineWidth', 2, 'DisplayName', ' Reproduction Number,R0');
plot(t, yd + 1.2, 'ko', 'LineWidth', 2, 'DisplayName', 'mean \delta');
plot(t_fit, yd_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Virulence\delta');
legend('show');
xlabel('days, t');
ylabel('mean values');
title('Relation betweeem R0 and \delta');
grid on;
% 
% plot(t, diff(yb), 'r-o', 'LineWidth', 2, 'DisplayName', 'mean \beta_s');
% grid on;
% hold on;
% plot(t, diff(yd), 'k--o', 'LineWidth', 2, 'DisplayName', 'mean \delta');
% legend(' transmission \beta_s',' virulence \delta')
% xlabel('days, t');
% ylabel('mean values');



%figure
% 
% plot(yb, yd, 'k--o', 'LineWidth', 2);
% xlabel('t');
% ylabel('mean values');
% grid on;
% Bar plot for bs_opt_values and del_opt_values separately
% Histograms for each wave for beta_s
figure;
for i = 1:6
    subplot(2, 3, i);
    histogram(bs_opt_values(:, i), 30);
    title(['Wave ' num2str(i) ' \beta_s']);
    xlabel('\beta_s values');
    ylabel('Frequency');
    grid on;
end
sgtitle('\beta_s Distribution for Each Wave');

% Histograms for each wave for delta
figure;
for i = 1:6
    subplot(2, 3, i);
    histogram(del_opt_values(:, i), 30);
    title(['Wave ' num2str(i) ' \delta']);
    xlabel('\delta values');
    ylabel('Frequency');
    grid on;
end
sgtitle('\delta Distribution for Each Wave');

%%
% Histograms for each wave for beta_s
figure;
for i = 1:7
    subplot(2, 4, i);
    histogram(bs_opt_values(:, i), 30, 'Normalization', 'probability');
    title(['Wave ' num2str(i) ' \beta_s']);
    xlabel('\beta_s values');
    ylabel('Probability');
    grid on;
end
sgtitle('\beta_s Distribution for Each Wave');

% Histograms for each wave for delta
figure;
for i = 1:7
    subplot(2, 4, i);
    histogram(del_opt_values(:, i), 30, 'Normalization', 'probability');
    title(['Wave ' num2str(i) ' \delta']);
    xlabel('\delta values');
    ylabel('Probability');
    grid on;
end
sgtitle('\delta Distribution for Each Wave');

% Histograms for each wave for xi
figure;
for i = 3:7
    subplot(2, 4, i-2);
    histogram(zi_opt_values(:, i), 30, 'Normalization', 'probability');
    title(['Wave ' num2str(i) ' \xi_{v}']);
    xlabel('\xi_{v} values');
    ylabel('Probability');
    grid on;
end
sgtitle('\xi_{v} Distribution for Each Wave');
%%
% Define time points
t = [92, 187, 365, 584, 765, 926, 1075]; 

% Estimates for bs
bs_estimates = bs_opt_values';
yb = mean(bs_estimates, 2)'; % Mean R0 values for each wave


% % Estimates for Rs
% Rs_estimates = Rs_opt_values';
% yr = mean(Rs_estimates, 2)'; % Mean R0 values for each wave

% Estimates for Rs
R0_estimates = R0_opt_values';
yr0 = mean(R0_estimates, 2)'; % Mean R0 values for each wave

% Estimates for delta
Delta_estimates = del_opt_values';
yd = mean(Delta_estimates, 2)'; % Mean delta values for each wave

% Estimates for zi
zi_estimates = zi_opt_values';
yz = mean(zi_estimates, 2)'; % Mean delta values for each wave


% Perform cubic spline interpolation
t_fit = linspace(min(t), max(t), 100);
yb_fit = spline(t, yb+1.5, t_fit);
% yr_fit = spline(t, yr, t_fit);
yr0_fit = spline(t, yr0, t_fit);
yd_fit = spline(t, yd+1.5, t_fit);
yz_fit = spline(t, yz+1.5, t_fit);


% % Plot the cubic spline fits
% figure;
% plot(t, yr, 'ro', 'LineWidth', 2, 'DisplayName', 'mean R0');
% hold on;
% plot(t_fit, yr_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Reproduction Number, R0');
% plot(t, yb+0.2, 'bo', 'LineWidth', 2, 'DisplayName', 'mean \beta_s ');
% plot(t_fit, yb_fit, 'b-', 'LineWidth', 2, 'DisplayName', 'transmission rate, \beta_s');
% plot(t, yd+0.4, 'ko', 'LineWidth', 2, 'DisplayName', 'mean \delta');
% plot(t_fit, yd_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Virulence \delta');
% plot(t, yz+0.3, 'go', 'LineWidth', 2, 'DisplayName', 'mean \xi_{v}');
% plot(t_fit, yz_fit, 'g--', 'LineWidth', 2, 'DisplayName', 'Vaccination \xi');


% Plot the cubic spline fits
figure;
plot(t, yr0, 'ro', 'LineWidth', 2, 'DisplayName', 'mean R0');
hold on;
plot(t_fit, yr0_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Reproduction Number, R0');
plot(t, yb+1.5, 'bo', 'LineWidth', 2, 'DisplayName', 'mean \beta_s ');
plot(t_fit, yb_fit, 'b-', 'LineWidth', 2, 'DisplayName', 'transmission rate, \beta_s');
plot(t, yd+1.5, 'ko', 'LineWidth', 2, 'DisplayName', 'mean \delta');
plot(t_fit, yd_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Virulence \delta');
plot(t, yz+1.5, 'go', 'LineWidth', 2, 'DisplayName', 'mean \xi_{v}');
plot(t_fit, yz_fit, 'g--', 'LineWidth', 2, 'DisplayName', 'Vaccination \xi');
legend('show');
xlabel('days, t');
ylabel('mean values');
title('Relation between R0,\beta_s and \delta');
grid on;
%%
% Histograms for each wave for beta_s
for i = 1:7
    figure;
    histogram(bs_opt_values(:, i), 30, 'Normalization', 'probability');
    title(['Wave ' num2str(i) ' \beta_s']);
    xlabel('\beta_s values');
    ylabel('Probability');
    grid on;
end

% Histograms for each wave for delta
for i = 1:7
    figure;
    histogram(del_opt_values(:, i), 30, 'Normalization', 'probability','FaceColor', [0.8500 0.3250 0.0980]);
    title(['Wave ' num2str(i) ' \delta']);
    xlabel('\delta values');
    ylabel('Probability');
    grid on;
end
%%
load('ba_opt_values.mat')
load('del_opt_values.mat')
% Estimates
ba1=ba_opt_values(:,1);
ba2=ba_opt_values(:,2);
ba3=ba_opt_values(:,3);
ba4=ba_opt_values(:,4);
ba5=ba_opt_values(:,5);
ba6=ba_opt_values(:,6);
ba7=ba_opt_values(:,7);


Betaa_estimate1 = ba1';% this is all simulates monte calro beta_s values for wave # 1
Betaa_estimate2 = ba2';
Betaa_estimate3 = ba3';
Betaa_estimate4 = ba4';
Betaa_estimate5 = ba5';
Betaa_estimate6 = ba6';
Betaa_estimate7 = ba7';

yba = [mean(Betaa_estimate1), mean(Betaa_estimate2),mean(Betaa_estimate3),mean(Betaa_estimate4),mean(Betaa_estimate5),mean(Betaa_estimate6),mean(Betaa_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
estimatesa = [Betaa_estimate1; Betaa_estimate2; Betaa_estimate3; Betaa_estimate4;Betaa_estimate5; Betaa_estimate6;Betaa_estimate7];
% Calculate the minimum and maximum bounds for shading
yba_min = min(estimatesa');
yba_max = max(estimatesa');

%%%%%%% estimates for delta

% Estimates
d1=del_opt_values(:,1);
d2=del_opt_values(:,2);
d3=del_opt_values(:,3);
d4=del_opt_values(:,4);
d5=del_opt_values(:,5);
d6=del_opt_values(:,6);
d7=del_opt_values(:,7);


Delta_estimate1 = d1';% this is all simulates monte calro delta values for wave # 1
Delta_estimate2 = d2';
Delta_estimate3 = d3';
Delta_estimate4 = d4';
Delta_estimate5 = d5';
Delta_estimate6 = d6';
Delta_estimate7 = d7';



yd = [mean(Delta_estimate1), mean(Delta_estimate2),mean(Delta_estimate3),mean(Delta_estimate4),mean(Delta_estimate5),mean(Delta_estimate6),mean(Delta_estimate7)]; % this is the mean beta_s values for each wave
% Combine estimates into a matrix
Destimates = [Delta_estimate1;Delta_estimate2;Delta_estimate3;Delta_estimate4;Delta_estimate5;Delta_estimate6;Delta_estimate7];

% Calculate the minimum and maximum bounds for shading
yd_min = min(Destimates');
yd_max = max(Destimates');

figure
sbetaa=diff(yba);
sdelta=diff(yd);
diffbad=[sbetaa' sdelta'];
bar(diffbad)
set(gca, 'XTickLabel', {'w1-w2', 'w2-w3', 'w3-w4', 'w4-w5', 'w5-w6','w6-w7'});
ylabel('Evolution in mean \beta_{a} and mean \delta values');
title('Relation between transmission and virulence');
legend('Change in \beta_a', 'Change in \delta', 'Location', 'Best');