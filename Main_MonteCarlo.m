clc
clear
load('COVIDSIRData.mat')
% Define initial parameters and their bounds
scale = 10^6;
N1= 331502651/scale;
global ta tb pp IC
 % Define the scale, adjust as needed

%a_min =23582486 * 1 / (67.7 * 350);a_max = 23582486 * 1 / (67.7 * 120);
a_min =9500/scale;%a_min/scale;
a_max =11670/scale;%a_max/scale;
w_min = 0.0004; w_max = 0.0100;
alpha_min = 0.00371; alpha_max =0.0125;
%bs_min =0.0031 ; bs_max = 0.6577;   %actual
bs_min =0.0031 ; bs_max = 0.3577;
%ba_min = 0.0021; ba_max = 0.5650;   %actual
ba_min = 0.00211; ba_max = 0.4650;
zi_min = 0.0024; zi_max = 0.0600;
e_min = 0.6; e_max = 0.94;
u_min = 0.3349*10^(-4); u_max =0.3915*10^(-4);
sigma_min = 0.06321; sigma_max =0.5208;  %actual
%sigma_min = 0.1532; sigma_max =0.5208;
r_min = 0.02341; r_max = 0.6534;   %actual
%r_min = 0.2341; r_max = 0.6534;
eta_min = 0.03; eta_max = 0.221;  %actual
%eta_min = 0.086; eta_max = 0.321;
%del_min = 0.0015; del_max = 0.5889;
%del_min = 0.00015; del_max = 0.001889; %actual
del_min = 0.00015; del_max = 0.05889;
phi_min = 0.02; phi_max = 0.156;   %actual
%phi_min = 0.046; phi_max = 0.256;

%%%%%%%%%%%%  here we set the day ta abd tb according to each bump. 

% 1• SARS COV-2 (original strain) : Day 1 to Day 180 (January 20, 2020 – June 18, 2020)
% 2• Alpha variant: Day 181 to Day 400 (June 18, 2020 – February 23, 2021)
% 3• Delta variant: Day 401 to Day 650 (February 23, 2021 – October 31, 2021)
% 4• Omicron variant: Day 651 to Day 850 (October 31, 2021 – May 19, 2022)
% 5• BA.2 variant: Day 851 to Day 1000 (May 19, 2022 – October 16, 2022)
% 6 XBB variant: Day 1001 to Day 1157 (October 16, 2022 – March 3, 2023, approximate)

numsim=100;
tVec=zeros(300,numsim);  % you need to increase 250 if you get some error due to large size of t vector
yVec=zeros(300,numsim);
casenum=1;
%%%%%%%%%%%%% 
for i=1:numsim

i
if casenum ==1
ta=30;
tb=155;

zi_min = 0; zi_max = 0;
w_min = 0; w_max = 0;
e_min = 0; e_max = 0;

pp=0.2; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
%IC = [N1- 1.3853, 0, 1200/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
IC = [30.998 , 0, 0.001, pp*I1data(1) , I1data(1), daily_recovered(ta)/scale];
elseif casenum ==11
ta=155;
tb=230;
zi_min = 0; zi_max = 0;
w_min = 0; w_max = 0;
e_min = 0; e_max = 0;
sigma_min = 0.4321; sigma_max = 1.4541;
r_min = 0.6341; r_max = 0.9534;
pp=0.3; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
%IC = [N1- 1.3853, 0, 1200/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
IC = [68.56, 0, 0.001, pp*I1data(1) , I1data(1), daily_recovered(ta)/scale];

elseif casenum ==2
ta=230;
tb=500;
pp=0.17; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
IC = [80.1174, 0.001/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
elseif casenum ==3
ta=500;
tb=665;
pp=0.11; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
IC = [110.1174, 0.001/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
elseif casenum ==4
ta=680;
tb=850;
pp=0.3; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
IC = [N1- 1.3853, 0.001/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
elseif casenum ==5
ta=851;
tb=1000;
%sigma_min = 0.8321; sigma_max = 2.4541;
%bs_min =0.009 ; bs_max = 0.3577;
%ba_min = 0.005; ba_max = 0.1650;
pp=0.4; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
IC = [50.112, 0.01/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
elseif casenum ==6
ta=1000;
tb=1150;
%sigma_min = 0.2321; sigma_max = 1.4541;
%r_min=0.0021; r_max=0.0551;
%bs_min =0.009 ; bs_max = 0.6577;
%ba_min = 0.005; ba_max = 0.5650;
pp=0.33; %percent aymptotic
I1data = daily_infectious(ta:tb,1)/scale;
IC = [30.561, 0.01/scale, 0.0002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];
end


% Random initial parameter values within bounds
a = a_min + (a_max - a_min) * rand;
w = w_min + (w_max - w_min) * rand;
alpha = alpha_min + (alpha_max - alpha_min) * rand;
bs = bs_min + (bs_max - bs_min) * rand;
ba = ba_min + (ba_max - ba_min) * rand;
zi = zi_min + (zi_max - zi_min) * rand;
e = e_min + (e_max - e_min) * rand;
u = u_min + (u_max - u_min) * rand;
sigma = sigma_min + (sigma_max - sigma_min) * rand;
r = r_min + (r_max - r_min) * rand;
eta = eta_min + (eta_max - eta_min) * rand;
del = del_min + (del_max - del_min) * rand;
phi = phi_min + (phi_max - phi_min) * rand;

%ta=500;
%tb=650;
%pp=0.4; %percent aymptotic
% Load actual infection data (replace 'A' with your actual data)
Tdata =1:tb-ta+1; %D1(ta:tb,1);
Tdata=Tdata';
%I1data = daily_infectious(ta:tb,1)/scale;

% Initial conditions
%IC = [N1 - 1.3853, 0.001/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];


% Time span
tspan =[Tdata(1) Tdata(end)];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-5 1e-4 1e-5 1e-4 1e-4]);

% Parameter bounds for optimization
LB =  [a_min, w_min, alpha_min, bs_min, ba_min, zi_min, e_min, u_min, sigma_min, r_min, eta_min, del_min, phi_min];
UB =  [a_max, w_max, alpha_max, bs_max, ba_max, zi_max, e_max, u_max, sigma_max, r_max, eta_max, del_max, phi_max];


% Define the optimization options
opts = optimset('MaxFunEvals', 800, 'MaxIter', 200);

% Initial guess for the parameters
par0 = [a, w, alpha, bs, ba, zi, e, u, sigma, r, eta, del, phi];

%%%% here we save the initial parameter values in a matrix. Each row is a
%%%% single simulation
par0vec(i,:)=par0;

% Perform the optimization
%[par_opt, fval] = fmincon(@(par) objFun(par, Tdata, I1data, IC, tspan, N), par0, [], [], [], [], LB, UB, [], opts)

[par_opt, fval] = fmincon(@objecFun, par0, [], [], [], [], LB, UB, [], opts);
% Extract optimized parameters
RSS=fval;
a_opt = par_opt(1);
w_opt = par_opt(2);
alpha_opt = par_opt(3);
bs_opt = par_opt(4);
ba_opt = par_opt(5);
zi_opt = par_opt(6);
e_opt = par_opt(7);
u_opt = par_opt(8);
sigma_opt = par_opt(9);
r_opt = par_opt(10);
eta_opt = par_opt(11);
del_opt = par_opt(12);
phi_opt = par_opt(13);

% Define the optimized parameters array
opt_par = [a_opt, w_opt, alpha_opt, bs_opt, ba_opt, zi_opt, e_opt, u_opt, sigma_opt, r_opt, eta_opt, del_opt, phi_opt];

%%%% here we save the estimated (optimized) parameter values in a matrix. Each row is asingle simulation
opt_parvec(i,:)=opt_par;

%%%%% please enter R0 exprestion here

R0=((sigma_opt*(u_opt+w_opt+zi_opt*(1-e_opt)))/((u_opt+sigma_opt)*(u_opt+w_opt+zi_opt)))...
    *(((r_opt*ba_opt)/(u_opt+eta_opt))+(((1-r_opt)*bs_opt)/(u_opt+del_opt+phi_opt)))

R0s=((sigma_opt*(u_opt+w_opt+zi_opt*(1-e_opt)))/((u_opt+sigma_opt)*(u_opt+w_opt+zi_opt)))...
    *(((1-r_opt)*bs_opt/(u_opt+del_opt+phi_opt)))

R0a=((sigma_opt*(u_opt+w_opt+zi_opt*(1-e_opt)))/((u_opt+sigma_opt)*(u_opt+w_opt+zi_opt)))...
    *((r_opt*ba_opt/(u_opt+eta_opt)))

%%%% here we save the estimated R0 values in a vector. 
R0vec(i)=R0;
R0avec(i)=R0a;
R0svec(i)=R0s;
%%%%% here we compure R^2 value

% R^2	=	coefficient of determination
% RSS	=	sum of squares of residuals
% TSS	=	total sum of squares

Iavg=mean(I1data);
TSS=sum((I1data-Iavg).^2);

Rsqure=1-(RSS/TSS)

goodnessOfFit(i,:)=[RSS TSS Rsqure];

% Run the ODE solver with optimized parameters
%[t_opt, y_opt] = ode45(@(t, p) sveair_model(t, p, opt_par, N), tspan, IC);
[t_opt, p_opt] = ode45(@sveair_model, tspan, IC,options,opt_par);

te=length(t_opt);
tVec(1:te,i)=t_opt;
yVec(1:te,i)=p_opt(:,5);


% % Plot results with optimized parameters
%figure;
%plot(t_opt, p_opt(:,5), '-r', 'LineWidth', 2.5); % Model infected people
%hold on;
%plot(Tdata, I1data, 'ob'); % Actual infected data
% xlabel('Time (days)');
% ylabel('Number of Infected Individuals per Million');
% legend('Model Solution', 'Actual Data');
% title('Model Solution Vs. Actual Infected Individual');
% grid on;

end

R0vec=R0vec';

R0avec=R0avec';
R0svec=R0svec';

Mestimate=[par0vec opt_parvec R0vec R0avec R0svec goodnessOfFit];
TYestimate=[tVec yVec];
writematrix(Mestimate,'Mcase11.xls')
%writematrix(TYestimate,'TYcase2.xls')
% this is the heading of your matrix. You should manually add it:

%a, w, alpha, bs, ba, zi, e, u, sigma, r, eta, del, phi, a_opt, 
%w_opt, alpha_opt, bs_opt, ba_opt, zi_opt, e_opt, u_opt, sigma_opt, r_opt, eta_opt, del_opt, phi_opt, R0, R0a R0s  RSS TSS Rsqure


stata_opt=datastats(Mestimate(:,14))
statw_opt=datastats(Mestimate(:,15)) 
statalpha_opt=datastats(Mestimate(:,16)) 
statbs_opt=datastats(Mestimate(:,17)) 
statba_opt=datastats(Mestimate(:,18))
statzi_opt=datastats(Mestimate(:,19)) 
state_opt=datastats(Mestimate(:,20))
statu_opt=datastats(Mestimate(:,21)) 
statsigma_opt=datastats(Mestimate(:,22)) 
statr_opt=datastats(Mestimate(:,23)) 
stateta_opt=datastats(Mestimate(:,24)) 
statdel_opt=datastats(Mestimate(:,25))
statphi_opt=datastats(Mestimate(:,26)) 

R0stat=datastats(R0vec)
R0astat=datastats(R0avec)
R0sstat=datastats(R0svec)
RSSstat=datastats(goodnessOfFit(:,1))
TSSstat=datastats(goodnessOfFit(:,2))
Rsqurestat=datastats(goodnessOfFit(:,3))
%  here we want to comppare the distributation of etiamted beta values
x1=Mestimate(:,14);
x2=Mestimate(:,15);
x3=Mestimate(:,16);
x4=Mestimate(:,17);
x5=Mestimate(:,18);
x6=Mestimate(:,19); 
x7=Mestimate(:,20);
x8=Mestimate(:,21);
x9=Mestimate(:,22);
x10=Mestimate(:,23); 
x11=Mestimate(:,24); 
x12=Mestimate(:,25);
x13=Mestimate(:,26);
figure(1)
h1=histogram(x1,'Normalization','probability', 'DisplayName','recruitment rate');
legend('show')
figure(2)
h2=histogram(x2,'Normalization','probability','DisplayName','vaccination wanning rate');
legend('show')
figure(3)
h3=histogram(x3,'Normalization','probability', 'DisplayName','reinfection rate');
legend('show')
figure(4)
h4=histogram(x4,'Normalization','probability','DisplayName','symptomatic transmission rate');
legend('show')
figure(5)
h5=histogram(x5,'Normalization','probability','DisplayName','asymptomatic transmission rate');
legend('show')
figure(6)
h6=histogram(x6,'Normalization','probability','DisplayName','Vaccination rate');
legend('show')
figure(7)
h7=histogram(x7,'Normalization','probability','DisplayName','Vaccine efficacy rate');
legend('show')
figure(8)
h8=histogram(x8,'Normalization','probability','DisplayName','Natural Mortality rate');
legend('show')
figure(9)
h9=histogram(x9,'Normalization','probability','DisplayName','Rate of exposed individual become infected');
legend('show')
figure(10)
h10=histogram(x10,'Normalization','probability','DisplayName','Ratio of asymptomatic infected from expose');
legend('show')
figure(11)
h11=histogram(x11,'Normalization','probability','DisplayName','recovered rate in asymptomatic class');
legend('show')
figure(12)
h12=histogram(x12,'Normalization','probability', 'DisplayName','infected mortality rate');
legend('show')
figure(13)
h13=histogram(x13,'Normalization','probability','DisplayName','recovered rate in symptomatic class');
legend('show')
% y = Mestimate(:,18);
% h1 = histogram(x);
% hold on
% h2 = histogram(y);
% h1.Normalization = 'probability';
% h1.BinWidth = 0.25;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.25;

%plot mean and standard deviation using Simbiology
% eeResults = sbioelementaryeffects(sveair_model,par0,yVec,NumberSamples=100,ShowWaitbar=true);
% plotData(eeResults,ShowMedian=true,ShowMean=false)
%sbioplot(tVec,yVec);


% Example vector
t_orig = tVec(:,1);

% Find the index of the last non-zero element
lastNonZeroIndex = find(t_orig, 1, 'last');

% Trim the vector to exclude zero entries at the bottom
trimmedt_orig = t_orig(1:lastNonZeroIndex);


for i=1:numsim

lastNonZeroIndexT = find(tVec(:,i), 1, 'last');
trimmedt_Vec = tVec(1:lastNonZeroIndexT,i);
lastNonZeroIndexY = find(yVec(:,i), 1, 'last');
trimmedy_Vec = yVec(1:lastNonZeroIndexY,i);
y_Vec_interp(:,i) = interp1(trimmedt_Vec, trimmedy_Vec, trimmedt_orig, 'linear', 'extrap');
end


y_min = min(y_Vec_interp,[],2);
y_max = max(y_Vec_interp,[],2);


figure(14);


% Plot the shaded region
fill([trimmedt_orig' fliplr(trimmedt_orig')], [y_min' fliplr(y_max')], [0.30 0.75 0.93], 'EdgeColor', 'none', 'DisplayName', 'Estimate Range','LineWidth', 1.5);
hold on

plot(Tdata, I1data, '.b'); % Actual infected data
xlabel('Time (days)');
ylabel('Number of Infected Individuals per Million');
grid on;
hold off;
legend('Estimate Rang','Covid-19 Data');

% figure
% 
% boxplot(R0vec)



