function [bite, piVH, piHV, k, p, t_peak, a, LambdaF, LambdaC, rF, rA, Error2] = OROV_MultiStart(data,SA0,IF0,SH0,IH0,SC0,IC0,IHC0,tspan,lb,ub,w,nu,NoStartPoints)

%% Upper and Lower bounds for parameters you are fitting.

%Parameter vector we are approximating:
% z = [b, piVH, piHV k, p, t_peak, a, LambdaF, LambdaC, rF, rA]

LowerBounds=lb;      %Lowerbounds for the parameters you are estimating
UpperBounds=ub;       %Upperbounds for the parameters you are estimating

xstart=.5*(LowerBounds+UpperBounds);                            %What initial parameter values you want to start with

%% MultiStart and fmincon - Fitting Part - Parallelization

% Here we set-up the optimization problem, specifying we will use fmincon
% as the local solver, and the what model we want to minimize along with
% the specific measure down below in the OROV_RUN_ODE15 function. We give
% initial conditions and the bounds.

SC0= (7*510/3)*(1/((426299/122.57)*(1.21+3.92)))*2054713;  %%initial number of susceptible city midges

%creating new initial values to include estimation of initial value of
%susceptible forest midges and infected sloths
new_initialvalues = @(rF, rA) [SA0-rA*SA0 rA*SA0 rF*SC0 IF0 SH0 IH0 SC0 IC0 IHC0]';

problem = createOptimProblem('fmincon','objective',@(z) OROV_RUN_ODE15(z,data,new_initialvalues(z(10),z(11)), tspan,w,nu)...
    ,'x0',xstart,'lb',LowerBounds,'ub',UpperBounds);%,'Aineq',A,'bineq',b)%,'Aeq',Aeq,'beq',beq);

% problem.options = optimoptions(problem.options,'TolFun',1e-16, 'TolX', 1e-16);%'MaxFunEvals',9999,'MaxIter',9999);%,'TolCon',0)
problem.options = optimoptions(problem.options,'MaxFunEvals',9999,'MaxIter',9999);%,'TolCon',0)
%problem.options = optimoptions(problem.options,'MaxFunEvals',inf,'MaxIter',inf,'TolFun',1e-10,'TolCon',0,'TolX',0,'MaxFunEvals',999999)

numstartpoints=NoStartPoints;                               % How many runs do you want?

% %  ms=MultiStart('Display','iter');                       %defines a multistart problem without parallel

ms=MultiStart('UseParallel',false,'Display','iter');         %defines a parallel multistart problem

%parpool %accesses the cores for parallel on your computer (laptop goes for 2-8, can be more specific)

[b,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  %runs the multistart

% the following takes solutions from manymins and makes a matrix out of them


for i=1:length(manymins)
    OROVParameters(i,:)=manymins(i).X;       %what are the parameter values
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;            %the minimization error
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;      %how "good" is the solution, we want 1 or 2.
end


%delete(gcp('nocreate'))  %turns off the parallel feature


% %% Plot the "best" solution
% infectiondata = [0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	1	0	1	1	6	10	32	53	124	178	227	229	321	357	311	270	197	254	197	175	119	111	59	72	61	52	41	41	24	17	15 5];
% % cumulativedata = [0	0	0	0	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	3	4	5	5	6	7	13	23	55	108	232	410	637	866	1187	1544	1855	2125	2322	2576	2773	2948	3067	3178	3237	3309	3370	3422	3463	3504	3528	3545	3560];
% 
% infectiondata_asymp = infectiondata/0.4;
% % 1% noise
% %noise = 0.01.*infectiondata_asymp.*randn(size(infectiondata_asymp));
% % 10% noise
% noise = 0.1.*infectiondata_asymp.*randn(size(infectiondata_asymp));
% 
% infectiondata_asymp = infectiondata_asymp + noise;
% 
%Outputs state variables for "best" fit
%New initial values
ini_val = [SA0-OROVParameters(1,11)*SA0 OROVParameters(1,11)*SA0 OROVParameters(1,10)*SC0 IF0 SH0 IH0 SC0 IC0 IHC0];

[t,y] = ode15s(@(t,y) OROV_Model(t,y,OROVParameters(1,:)),tspan,ini_val);
    SA = y(:,1);
    IA = y(:,2);
    SF = y(:,3);
    IF = y(:,4);
    SH = y(:,5);
    IH = y(:,6);
    SC = y(:,7);
    IC = y(:,8);
    IHC = y(:,9);
% 
% 
%     plot(tspan,IH)   %plots "best fit"
%     hold on
% 
% %plot the actual data
% scatter(tspan,infectiondata_asymp,'filled')
% legend('Infected Humans', 'Location','northwest')
% xlabel('Weeks')
% ylabel('Infections')
% title('Best Fit (continuous curve) with Weekly Data (dots)')


% plot(tspan,IHC)   %plots "best fit"
%     hold on
% 
% %plot the actual data
% scatter(tspan,cumulativedata)
% legend('Infected Humans', 'Location','northwest')
% xlabel('Weeks')
% ylabel('Infections')
% title('Best Fit (continuous curve) with Cumulative Data (dots)')

%Error for weekly data
Error2 = rl2Err(IH, data); % Calculate the l2 Error

% %Error for cumulative data
% Error2 = rl2Err(IHC, data); % Calculate the l2 Error

%%define OROVParameters
bite = OROVParameters(1,1);
piVH= OROVParameters(1,2);
piHV = OROVParameters(1,3);
k = OROVParameters(1,4);
p = OROVParameters(1,5);
t_peak = OROVParameters(1,6);
a = OROVParameters(1,7);
LambdaF = OROVParameters(1,8);
LambdaC = OROVParameters(1,9);
rF = OROVParameters(1,10);
rA = OROVParameters(1,11);
