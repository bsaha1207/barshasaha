clear
clc
%%%%%%this code apply interpolation CHIP to the histograms of beta and delta

%t = [1, 2, 3, 4, 5, 6, 7];  %this is the epidemic wave number
% t = [92, 187, 365, 584, 765, 926, 1075]; 
t = [187, 365, 584, 765, 926, 1075]; 
%y = [2, 4, 9, 16]; % this is the mean beta_s values for each wave
load('bs_opt_values.mat')
load('ba_opt_values.mat')
load('del_opt_values.mat')
load('zi_opt_values.mat')
load('COVIDSIRData.mat')

% Estimates
b1=bs_opt_values(:,1);
b2=bs_opt_values(:,2);
b3=bs_opt_values(:,3);
b4=bs_opt_values(:,4);
b5=bs_opt_values(:,5);
b6=bs_opt_values(:,6);
b7=bs_opt_values(:,7);

% Estimates of ba
ba1=ba_opt_values(:,1);
ba2=ba_opt_values(:,2);
ba3=ba_opt_values(:,3);
ba4=ba_opt_values(:,4);
ba5=ba_opt_values(:,5);
ba6=ba_opt_values(:,6);
ba7=ba_opt_values(:,7);


%%%%%%% estimates for delta

% Estimates
d1=del_opt_values(:,1);
d2=del_opt_values(:,2);
d3=del_opt_values(:,3);
d4=del_opt_values(:,4);
d5=del_opt_values(:,5);
d6=del_opt_values(:,6);
d7=del_opt_values(:,7);

% Estimates
z1=zi_opt_values(:,1);
z2=zi_opt_values(:,2);
z3=zi_opt_values(:,3);
z4=zi_opt_values(:,4);
z5=zi_opt_values(:,5);
z6=zi_opt_values(:,6);
z7=zi_opt_values(:,7);


%[xFine1, fPchip1]= InterpolDistfitter(b1)

[xb1, fb1]=InterpolDistfitter(b1);
[xb2, fb2]=InterpolDistfitter(b2);
[xb3, fb3]=InterpolDistfitter(b3);
[xb4, fb4]=InterpolDistfitter(b4);
[xb5, fb5]=InterpolDistfitter(b5);
[xb6, fb6]=InterpolDistfitter(b6);
[xb7, fb7]=InterpolDistfitter(b7);

[xba1, fba1]=InterpolDistfitter(ba1);
[xba2, fba2]=InterpolDistfitter(ba2);
[xba3, fba3]=InterpolDistfitter(ba3);
[xba4, fba4]=InterpolDistfitter(ba4);
[xba5, fba5]=InterpolDistfitter(ba5);
[xba6, fba6]=InterpolDistfitter(ba6);
[xba7, fba7]=InterpolDistfitter(ba7);


[xd1, fd1]=InterpolDistfitter(d1);
[xd2, fd2]=InterpolDistfitter(d2);
[xd3, fd3]=InterpolDistfitter(d3);
[xd4, fd4]=InterpolDistfitter(d4);
[xd5, fd5]=InterpolDistfitter(d5);
[xd6, fd6]=InterpolDistfitter(d6);
[xd7, fd7]=InterpolDistfitter(d7);

[xz3, fz3]=InterpolDistfitter(z3);
[xz4, fz4]=InterpolDistfitter(z4);
[xz5, fz5]=InterpolDistfitter(z5);
[xz6, fz6]=InterpolDistfitter(z6);
[xz7, fz7]=InterpolDistfitter(z7);



%%%%%%%%%%%%%%%%%%%%%%%%%%%% for beta
% Define line styles and colors
lineStyles = {'-', '--', ':', '-.', '-', '--', ':'};
%lineStyles = { '--', ':', '-.', '-', '--', ':'};
colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
%colors = {'r', 'g', 'b', 'c', 'm', 'y'};

% Combine line styles and colors
lineSpec = strcat(colors, lineStyles);

% Create a new figure
figure;

% Plot all time series with different line styles and colors
 % plot(xb2,fb2, lineSpec{2}, ...
 %     xb3,fb3, lineSpec{3}, ...
 %     xb4,fb4, lineSpec{4}, ...
 %     xb5,fb5, lineSpec{5}, ...
 %     xb6,fb6, lineSpec{6}, ...
 %     xb7,fb7, lineSpec{7}, ...
 %     'LineWidth', 1.5);
 % 

 plot(xb1,fb1, lineSpec{1}, ...
     xb2,fb2, lineSpec{2}, ...
     xb3,fb3, lineSpec{3}, ...
     xb4,fb4, lineSpec{4}, ...
     xb5,fb5, lineSpec{5}, ...
     xb6,fb6, lineSpec{6}, ...
     xb7,fb7, lineSpec{7}, ...
     'LineWidth', 1.5);

% Add legend
legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
%legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6'});

% Display grid
grid on;


% Add labels and title
xlabel('Symptomatic Transmission rates (\beta_s) ');
ylabel('Probability Density Estimations of \beta_s');
%title('Kernel Density Estimation');
%legend('Data Histogram', 'Kernel Density Estimate');

figure;
plot(xba1,fba1, lineSpec{1}, ...
     xba2,fba2, lineSpec{2}, ...
     xba3,fba3, lineSpec{3}, ...
     xba4,fba4, lineSpec{4}, ...
     xba5,fba5, lineSpec{5}, ...
     xba6,fba6, lineSpec{6}, ...
     xba7,fba7, lineSpec{7}, ...
     'LineWidth', 1.5);

% Add legend
legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
%legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6'});

% Display grid
grid on;


% Add labels and title
xlabel('Asymptomatic Transmission rates (\beta_a) ');
ylabel('Probability Density Estimations of \beta_a');
%title('Kernel Density Estimation');
%legend('Data Histogram', 'Kernel Density Estimate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%% for delta
% Create a new figure
figure;

% Plot all time series with different line styles and colors
% plot(xd1,fd1, lineSpec{1}, ...
%      xd2,fd2/(9*10^2), lineSpec{2}, ...
%      xd3,fd3, lineSpec{3}, ...
%      xd4,fd4, lineSpec{4}, ...
%      xd5,fd5, lineSpec{5}, ...
%      xd6,fd6, lineSpec{6}, ...
%      xd7,fd7, lineSpec{7}, ...
%      'LineWidth', 1.5);

plot(xd1,fd1, lineSpec{1}, ...
    xd2,fd2, lineSpec{2}, ...
     xd3,fd3, lineSpec{3}, ...
     xd4,fd4, lineSpec{4}, ...
     xd5,fd5, lineSpec{5}, ...
     xd6,fd6, lineSpec{6}, ...
     xd7,fd7, lineSpec{7}, ...
     'LineWidth', 1.5);


% Add legend
legend({'w1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
%legend({'w1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6'});

% Display grid
grid on;
xlim([0 0.55]);


% Add labels and title
xlabel('Virulance rates (\delta) ');
ylabel('Probability Density Estimations of \delta');
%title('Kernel Density Estimation');
%legend('Data Histogram', 'Kernel Density Estimate');

% Create a new figure
figure;

% Plot all time series with different line styles and colors
% plot(xd1,fd1, lineSpec{1}, ...
%      xd2,fd2/(9*10^2), lineSpec{2}, ...
%      xd3,fd3, lineSpec{3}, ...
%      xd4,fd4, lineSpec{4}, ...
%      xd5,fd5, lineSpec{5}, ...
%      xd6,fd6, lineSpec{6}, ...
%      xd7,fd7, lineSpec{7}, ...
%      'LineWidth', 1.5);

plot(xd3,fd3, lineSpec{3}, ...
     xd4,fd4, lineSpec{4}, ...
     xd5,fd5, lineSpec{5}, ...
     xd6,fd6, lineSpec{6}, ...
     xd7,fd7, lineSpec{7}, ...
     'LineWidth', 1.5);


% Add legend
legend({'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
%legend({'w1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6'});

% Display grid
grid on;
xlim([0 0.5]);

% Add labels and title
xlabel('Vacciantion rates (\xi_{v}) ');
ylabel('Probability Density Estimations of \xi_{v}');


%%
% Load `ps` and `pa` values for each wave from separate files
ps1 = readmatrix('ps_values_1.xlsx');
ps2 = readmatrix('ps_values_2.xlsx');
ps3 = readmatrix('ps_values_3.xlsx');
ps4 = readmatrix('ps_values_4.xlsx');
ps5 = readmatrix('ps_values_5.xlsx');
ps6 = readmatrix('ps_values_6.xlsx');
ps7 = readmatrix('ps_values_7.xlsx');

pa1 = readmatrix('pa_values_1.xlsx');
pa2 = readmatrix('pa_values_2.xlsx');
pa3 = readmatrix('pa_values_3.xlsx');
pa4 = readmatrix('pa_values_4.xlsx');
pa5 = readmatrix('pa_values_5.xlsx');
pa6 = readmatrix('pa_values_6.xlsx');
pa7 = readmatrix('pa_values_7.xlsx');

% Define line styles and colors
lineStyles = {'-', '--', ':', '-.', '-', '--', ':'};
colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
lineSpec = strcat(colors, lineStyles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% for `ps` - Probability of Symptomatic Infection

[xps1, fps1] = InterpolDistfitter(ps1);
[xps2, fps2] = InterpolDistfitter(ps2);
[xps3, fps3] = InterpolDistfitter(ps3);
[xps4, fps4] = InterpolDistfitter(ps4);
[xps5, fps5] = InterpolDistfitter(ps5);
[xps6, fps6] = InterpolDistfitter(ps6);
[xps7, fps7] = InterpolDistfitter(ps7);

figure;
plot(xps1, fps1, lineSpec{1},...
xps2, fps2, lineSpec{2},....
xps3, fps3, lineSpec{3}, ....
xps4, fps4, lineSpec{4},...
xps5, fps5, lineSpec{5},....
xps6, fps6, lineSpec{6},....
xps7, fps7, lineSpec{7},...
'Linewidth',1.5)

xlabel('Probability of Symptomatic Infection (p_s)', 'FontSize', 12);
ylabel('Probability Density Estimation of p_s', 'FontSize', 12);
legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
grid on;
title('Probability Density Estimation of p_s for Each Wave', 'FontSize', 14);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% for `pa` - Probability of Asymptomatic Infection

[xpa1, fpa1] = InterpolDistfitter(pa1);
[xpa2, fpa2] = InterpolDistfitter(pa2);
[xpa3, fpa3] = InterpolDistfitter(pa3);
[xpa4, fpa4] = InterpolDistfitter(pa4);
[xpa5, fpa5] = InterpolDistfitter(pa5);
[xpa6, fpa6] = InterpolDistfitter(pa6);
[xpa7, fpa7] = InterpolDistfitter(pa7);

figure;
plot(xpa1, fpa1, lineSpec{1},...
xpa2, fpa2, lineSpec{2},....
xpa3, fpa3, lineSpec{3}, ....
xpa4, fpa4, lineSpec{4},...
xpa5, fpa5, lineSpec{5},....
xpa6, fpa6, lineSpec{6},....
xpa7, fpa7, lineSpec{7},...
'Linewidth',1.5)
xlim([0 0.075])
xlabel('Probability of Asymptomatic Infection (p_a)', 'FontSize', 12);
ylabel('Probability Density Estimation of p_a', 'FontSize', 12);
legend({'w 1', 'w 2', 'w 3', 'w 4', 'w 5', 'w 6', 'w 7'});
grid on;
title('Probability Density Estimation of p_a for Each Wave', 'FontSize', 14);
hold off;