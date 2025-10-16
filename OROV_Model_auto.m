function dydt = OROV_Model_auto(t,y,z,fT)

% z = [b, piVH, piHV k, p, t_peak, a, LambdaF, LambdaC, rF, rA] 
bite=z(1);
piVH=z(2);
piHV=z(3);
k=z(4);
p=z(5);
t_peak=z(6);
a=z(7);
LambdaF=z(8);
LambdaC=z(9);

% Parameters from literature
muA = 7*1/9125;    % 25 year sloth life expectancy
muH = 7*1/26801.95;% 73.43 years amazonas (human) life expectancy (2022 BR census)
muV = 7*1/25;      % 25 day midge lifespan
gamma = 7*1/4.5; % 4.5 day recovery
Nc = (7*510/3)*(1/((426299/122.57)*(1.21+3.92)))*2054713;
Nf = Nc*z(10);

% Initiate DE variables
dydt = zeros(9,1);

% Differential Equations -----------------------
%y=reshape(y,[],8);
SA = y(1);
IA = y(2);
SF = y(3);
IF = y(4);
SH = y(5);
IH = y(6);
SC = y(7);
IC = y(8);
IHC = y(9);

%Precipdata
precipData = [0.4528643113	0.5388131874	0.3626029337	0.1937340159	0.2988101152	0.2576278929	0.2260352204	0.242688757	0.2166257792	0.4208245113	0.3000935577	0.1654524679	0.08193460933	0.1691139284	0.1724167435	0.4086885127	0.3930246753	0.2359257	0.1047895424	0.3148334961	0.420935451	0.4035179125	0.3813075079	0.3985096796	0.2340377033	0.523321374	0.3779442822	0.5786127053	0.6197574238	0.7686275452	0.8618833974	0.7067981804	0.6158471349	0.5049413097	0.6107744492	0.8920679885	0.7200471696	0.7631654475	0.4932924129	0.6747123427	0.7037940491	0.5790355789	0.7694555509	0.8273047842	0.735913573	0.5503221744	0.9029677053	0.7086093727	0.5391208542	0.5537307641	0.571046795	0.6158258004	0.427956544];

tspanP = 1:1:length(precipData);

%Seasonality Functions
exponential = exp(-((t-t_peak)/a)^2);
cosine = 1/2*((1+a) + (1-a)*cos(2*pi*(t-t_peak)/tspanP(end)));

denomAlp = (a-1)/(a+t_peak-2);
tScaled = t/tspanP(end);
betadist = (tScaled/denomAlp)^(a-1)*((1-tScaled)/(1-denomAlp))^(t_peak-1);

preciplinear = interp1(tspanP,precipData,t,'linear','extrap');

switch fT(3)
    case 0
        seasonForce = 1;
    case 1
        seasonForce = exponential;
    case 2
        seasonForce = cosine;
    case 3
        seasonForce = betadist;
    case 4
        seasonForce = preciplinear;
end

if fT(1) == 1
    betaFA = bite*piVH*(k*(SA+IA))/(k*(SA+IA)+(1-p)*(SH+IH))*seasonForce;
    betaAF = bite*piHV*(k*(SA+IA))/(k*(SA+IA)+(1-p)*(SH+IH))*seasonForce;
    betaFH = bite*piVH*((1-p)*(SH+IH))/(k*(SA+IA)+(1-p)*(SH+IH))*seasonForce;
    betaHF = bite*piHV*((1-p)*(SH+IH))/(k*(SA+IA)+(1-p)*(SH+IH))*seasonForce;
    betaCH = bite*piVH*seasonForce;
    betaHC = bite*piHV*seasonForce;
else
    betaFA = bite*piVH*(k*(SA+IA))/(k*(SA+IA)+(1-p)*(SH+IH))*1;
    betaAF = bite*piHV*(k*(SA+IA))/(k*(SA+IA)+(1-p)*(SH+IH))*1;
    betaFH = bite*piVH*((1-p)*(SH+IH))/(k*(SA+IA)+(1-p)*(SH+IH))*1;
    betaHF = bite*piHV*((1-p)*(SH+IH))/(k*(SA+IA)+(1-p)*(SH+IH))*1;
    betaCH = bite*piVH*1;
    betaHC = bite*piHV*1;
end

if fT(2) == 1
    ReprF = LambdaF*seasonForce;
    ReprC = LambdaC*seasonForce;
else
    ReprF = muV*Nf;
    ReprC = muV*Nc;
end

dydt(1) = muA*(SA+IA) - betaFA*IF*SA/(SA+IA) + gamma*IA - muA*SA;
dydt(2) = betaFA*IF*SA/(IA+SA) - gamma*IA - muA*IA;
dydt(3) = ReprF - (betaAF*IA/(SA+IA)+betaHF*IH/(IH+SH))*SF - muV*SF;
dydt(4) = (betaAF*IA/(IA+SA)+betaHF*IH/(IH+SH))*SF - muV*IF;
dydt(5) = muH*(IH+SH) - (betaFH*IF+betaCH*IC)*SH/(SH+IH) + gamma*IH - muH*SH;
dydt(6) = (betaFH*IF+betaCH*IC)*SH/(SH+IH) - gamma*IH - muH*IH;
dydt(7) = ReprC - betaHC*IH/(IH+SH)*SC - muV*SC;
dydt(8) = betaHC*IH/(SH+IH)*SC - muV*IC;
dydt(9) = (betaFH*IF+betaCH*IC)*SH/(SH+IH);

end

