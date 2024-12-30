function err = objecFun(par)
%UNTITLED Summary of this function goes here
scale = 10^6;
N1 = 331502651/scale;
%   Detailed explanation goes here
    a = par(1);
    w = par(2);
    alpha = par(3);
    bs = par(4);
    ba = par(5);
    zi = par(6);
    e = par(7);
    u = par(8);
    sigma = par(9);
    r = par(10);
    eta = par(11);
    del = par(12);
    phi = par(13);
    pa=[a,w,alpha,bs,ba,zi,e,u,sigma,r,eta,del,phi];
load('COVIDSIRData.mat')
global ta tb pp IC

% Load actual infection data (replace 'A' with your actual data)

Tdata =1:tb-ta+1; %D1(ta:tb,1);
Tdata=Tdata';
I1data = daily_infectious(ta:tb,1)/scale;
Cdata = [Tdata, I1data];
% Initial conditions
%IC = [N1 - 1.3853, 0.001/scale, 0.002/scale, pp*I1data(1), I1data(1), daily_recovered(ta)/scale];

% Time span
tspan =[Tdata(1) Tdata(end)];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-5 1e-4 1e-5 1e-4 1e-4]);

[t,p] =ode45(@sveair_model,Tdata',IC,options,pa);

value = (p(:,5)-Cdata(:,2)).^2; %+1*(y(:,2)-Cdata(:,3)).^2;
err= sum(value);
%I_model = interp1(t, y(:,5), Tdata);
    %err = sum((I_model - I1data).^2)

end
