% A Novel Mathematical Model of Seasonality in Oropouche Virus Transmission in Amazonas, Brazil
% Data fitting code
%%
clear
clc
close all
format long
%% Infection data for May 28th, 2023 - June 1st, 2024

% Regular Infection & Cummulative Data
infectiondata = [0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	1	0	1	1	6	10	32	53	124	178	227	229	321	357	311	270	197	254	197	175	119	111	59	72	61	52	41	41	24	17	15 5];
cumulativedata = [0	0	0	0	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	3	4	5	5	6	7	13	23	55	108	232	410	637	866	1187	1544	1855	2125	2322	2576	2773	2948	3067	3178	3237	3309	3370	3422	3463	3504	3528	3545	3560];

% Modified Infection Data (divide by 0.4 for asymmptomatic people)
infectiondata_asymp = infectiondata/0.4;

% Add noise if desired
% % 1% noise
% noise = 0.01.*infectiondata_asymp.*randn(size(infectiondata_asymp));
% 10% noise
% noise = 0.1.*infectiondata_asymp.*randn(size(infectiondata_asymp));
% infectiondata_asymp = infectiondata_asymp + noise;

%Time span
tspan = 1:1:length(infectiondata_asymp);       %%vector - how many timepoints you are simulating, will be inputed into ode15s

%Initial Values
IH0=infectiondata_asymp(1); %%initial number of infected humans
SH0=3.941e6-IH0;   %%initial number of susceptible humans
IF0=0;  %%initial number of infected forest midges
SC0= (7*510/3)*(1/((426299/122.57)*(1.21+3.92)))*2054713; %%initial number of susceptible city midges
IC0=0;  %%initial number of infected city midges
SA0=(1559256.365-11401)*1.7;  %%initial number of susceptible sloths
IHC0=cumulativedata(1); %%initial number of cumulative cases

initialvalues = [SA0 IF0 SH0 IH0 SC0 IC0 IHC0]'; %%intial value vector

%size must be sxn where s is the number of solution componants and 
%n is the number of initial conditions being solved for. Each column in the
%matrix then respresents one complete set of initial conditions.

%% Data-fitting

% Running the multistart algorithm

strParams = ["$bite$","$\pi_VH$","$\pi_HV","$\k$","$p$","$t_peak$","$a$", "$\Lambda_F$", "$\Lambda_C$","$rF$", "$rA$"];

% Lower and upper bound for parameters being fit
% [bite, pi_VH, pi_HV, k, p, t_peak, a, LambdaF, LambdaC, rF, rA]
lb=[0,   0, 0, 0, 0,  0, 0, 0, 0, 0, 0];
ub=[20, 1, 1, 20, 1,0, 0, 0, 0, 20, 1];

% Regularization weight
nu = 0.1;

n=2000;    % n trials of the multistart algorithm
runs = 20;  % number of times fmincon runs
paramsAndErr = zeros(n,length(strParams)+2);  % matrix holding values of each parameter and relative l2 error percent for each trial of the multistart alogrithm

for i=0:n-1
    i+1
    w = 10^randi([-3,4]); % Data fitting weight
    [paramsAndErr(i+1,1), paramsAndErr(i+1,2), paramsAndErr(i+1,3), paramsAndErr(i+1,4), paramsAndErr(i+1,5), paramsAndErr(i+1,6), paramsAndErr(i+1,7), paramsAndErr(i+1,8), paramsAndErr(i+1,9), paramsAndErr(i+1,10), paramsAndErr(i+1,11), paramsAndErr(i+1,12)] = OROV_MultiStart(infectiondata_asymp, SA0, IF0, SH0, IH0, SC0, IC0, IHC0, tspan, lb, ub,w,nu,runs);
    paramsAndErr(i+1,13) = w;
end

[min, minRow] = min(paramsAndErr(:,12)); % find run with the smallest error
minRow; % outputs the row number with the smallest error

% the "optimal" parameter values
bite = paramsAndErr(minRow,1)
pi_VH = paramsAndErr(minRow,2)
pi_HV = paramsAndErr(minRow,3)
k = paramsAndErr(minRow,4)
p = paramsAndErr(minRow,5)
t_peak = paramsAndErr(minRow,6)
a = paramsAndErr(minRow,7)
Lambda_F = paramsAndErr(minRow,8)
Lambda_C = paramsAndErr(minRow,9)
rF = paramsAndErr(minRow,10)
rA = paramsAndErr(minRow,11)
w = paramsAndErr(minRow,13)


writematrix(paramsAndErr, 'Parameters.csv')      % creates a .csv file of the parameter and relative l2 error percent values
params = [paramsAndErr(minRow,1), paramsAndErr(minRow,2), paramsAndErr(minRow,3), paramsAndErr(minRow,4), paramsAndErr(minRow,5), paramsAndErr(minRow,6), paramsAndErr(minRow,7), paramsAndErr(minRow, 8), paramsAndErr(minRow, 9), paramsAndErr(minRow, 10), paramsAndErr(minRow, 11)];
rA = params(11); rF = params(10);
initialValues = [initialvalues(1)-rA*initialvalues(1) rA*initialvalues(1) rF*initialvalues(5) initialvalues(2) initialvalues(3) initialvalues(4) initialvalues(5) initialvalues(6) initialvalues(7)];

%% Section 5: plotting
colors = 1/255*[0 0 255; 255 0 0];%0 150 104; 239 108 0];

fineMesh = linspace(1, tspan(end),1001); % for pretty pictures
[t,y]=ode15s(@(t,y) OROV_Model(t,y,params), fineMesh, initialValues);

%infected human pop (I_h)
figure(3)
plot(fineMesh,y(:,6),'LineWidth',1.5);
xlim([0 fineMesh(end)+1])
hold on
scatter(tspan, infectiondata_asymp,'filled','r')
xlabel('Time (weeks)')
ylabel('Infections')
legend('simulated', 'observed',Location='northeast')
hold off

%infected sloth pop (I_a)
figure(4)
plot(fineMesh, y(:,2), 'color', colors(1,:), 'LineWidth', 1.5)
legend('Infected Sloths',Location='northeast')
xlabel('Time (weeks)')
ylabel('Sloths')
xlim([0 fineMesh(end)+1])
% set(gca, 'yscale','log')

%infected forest midge pop (I_f)
figure(5)
plot(fineMesh, y(:,4),'color', colors(2,:), 'LineWidth',1.5) % fineMesh,y(:,7),'b--',fineMesh,y(:,8),'r--',
legend('Infected Forest Midge', Location='northeast') % ,'Susceptible City Midge','Infected City Midge'
xlabel('Time (weeks)')
ylabel('Midges')
xlim([0 fineMesh(end)+1])
% set(gca, 'yscale','log')

%infected city midge pop (I_c)
figure(6)
plot(fineMesh,y(:,8),'--','color', colors(2,:), 'LineWidth',1.5)
legend('Infected City Midge')
xlabel('Time (weeks)')
ylabel('Midges')
xlim([0 fineMesh(end)+1])

%all infected pops to show peak times
figure(7)
plot(fineMesh,y(:,6),'LineWidth',1.5)
hold on
plot(fineMesh,y(:,2),'color', colors(1,:),'LineWidth',1.5)
hold on
plot(fineMesh,y(:,4),'color', colors(2,:),'LineWidth',1.5)
hold on
plot(fineMesh,y(:,8),'--','color', colors(2,:),'LineWidth',1.5)
legend('Infected Humans','Infected Sloths','Infected Forest Midges', 'Infected City Midges')%,Location='northwest'
xlabel('Time (weeks)')
ylabel('Count')
xlim([0 fineMesh(end)+1])
% set(gca, 'yscale','log')
hold off

%logscale midge pops
figure(8)
plot(fineMesh, y(:,4),'color', colors(2,:), 'LineWidth',1.5) % fineMesh,y(:,7),'b--',fineMesh,y(:,8),'r--',
hold on
plot(fineMesh,y(:,8),'--','color', colors(2,:), 'LineWidth',1.5)
legend('Infected Forest Midge','Infected City Midge')
xlabel('Time (weeks)')
ylabel('Midges')
xlim([0 fineMesh(end)+1])
set(gca, 'yscale','log')
hold off

%cumulative infected human pop (I_hc)
figure(9)
plot(fineMesh,y(:,9),'LineWidth',1.5);
xlim([0 fineMesh(end)+1])
hold on
scatter(tspan, cumulativedata,'filled','r')
xlabel('Time (weeks)')
ylabel('Infections')
legend('simulated', 'observed',Location='northeast')
hold off