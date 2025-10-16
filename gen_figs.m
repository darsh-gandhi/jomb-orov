clear; close all; format long;

%% raw data -> matrix

infectiondata = [0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	1	0	1	1	6	10	32	53	124	178	227	229	321	357	311	270	197	254	197	175	119	111	59	72	61	52	41	41	24	17	15];
cumulativedata = [0	0	0	0	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	3	4	5	5	6	7	13	23	55	108	232	410	637	866	1187	1544	1855	2125	2322	2576	2773	2948	3067	3178	3237	3309	3370	3422	3463	3504	3528	3545	3560];

infectiondata=infectiondata/0.4;
cumulativedata=cumulativedata/0.4;

%line 12 = raw precip data
precipData=[4.849141905	6.867277143	4.443242857	2.224508571	2.93492	3.033474286	2.895797143	2.960337143	2.878477143	4.106557143	3.594805714	1.900137143	1.768621429	2.292584286	2.821328571	5.046128571	4.730245714	3.300805714	1.819934286	3.552162857	4.486657143	4.331888571	4.662574286	5.71312	3.168265714	5.83674	4.732437143	8.666571429	8.164225714	8.6239	9.152665714	8.595968571	8.48332	7.181268571	7.237161905	9.87606381	8.700545714	10.53928857	6.766277143	8.363442857	8.900294286	6.962757143	9.624942857	9.398892381	8.95343619	6.957922857	12.37784571	7.783937143	7.3135	8.214531429	6.572394286	8.00038];
precipData=precipData./max(precipData); % normalize data
precipData_3dayavg = [0.5411453173	0.4514268108	0.3650500456	0.2850490216	0.2500573413	0.2608244095	0.2421172901	0.2284499189	0.2933796825	0.312514616	0.2954568456	0.182493545	0.1388336685	0.1411550938	0.2500730615	0.3247099772	0.3458796293	0.2445799726	0.2185162462	0.2801861632	0.3797622865	0.4019202905	0.3944450333	0.3379516303	0.3852895856	0.3784344532	0.4932927872	0.5254381371	0.6556658914	0.7500894555	0.779103041	0.7281762376	0.6091955417	0.5771876313	0.6692612492	0.7409632024	0.7917602019	0.65883501	0.643723401	0.6239329349	0.6525139902	0.6840950596	0.7252653047	0.7775579694	0.7045135105	0.7297344842	0.7206330841	0.7168993107	0.600486997	0.5546328044	0.5802011198	0.5382763798	0.4912677752];


%normalized precip data
%precipData = [0.391759763	0.554803905	0.358967381	0.179716941	0.237110727	0.245072879	0.233950011	0.239164166	0.232550737	0.331766709	0.290422566	0.153511135	0.142886046	0.185216744	0.227933732	0.407674218	0.382154199	0.266670452	0.147031586	0.286977471	0.362474799	0.349971124	0.376687058	0.461560124	0.255962612	0.471547322	0.382331243	0.700167996	0.659583736	0.696720593	0.739439312	0.694464026	0.685363204	0.580171117	0.584686711	0.797882284	0.702912762	0.851463883	0.546644166	0.675678389	0.719050349	0.562517687	0.777594347	0.759331842	0.723343657	0.562127128	1	0.628860411	0.590854028	0.663647909	0.530980466	0.646346722];

%Time span
tspan = 1:1:length(infectiondata);

colors = 1/255*[0 0 255; 255 0 0; 0 150 255; 0 0.4470*255 0.7410*255];% 239 108 0]; % blue, red, light blue, matlab blue

xr_colors = 1/255 * [115 147 179; 137 207 240];

colororder({'k'})

figure(1)
yyaxis left
scatter(tspan,infectiondata,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))
xlabel('Time (weeks)')
xlim([0.5 52])
ylim([-5 925])
ylabel('Human Infections')

yyaxis right
plot(tspan,precipData,'LineWidth',1.5,'Color',colors(4,:))
ylabel('Normalized Precipitation','Rotation',270)
ylim([0 1])
% ylim([0 1.008333333333333])


xr = xregion(27,52,FaceColor=xr_colors(2,:),EdgeColor=xr_colors(1,:));
% xr.FaceAlpha=0.5;
xr.EdgeAlpha=0.5;
xr.LineWidth=2;
legend('OROV Incidence','Precipitation','Rainy Season','Location', 'northwest')

%title('OROV Human Infection Data June 2023 - June 2024')
set(gca,'fontsize',12)

%% parameters

% params = [(1) bite, (2) pi_VH, (3) pi_HV, (4) k, (5) p, (6) t_peak, (7) a, (8) LambdaF, (9) LambdaC, (10) rF, (11) rA];
% %no seasonality
x=[7.019493013	0.238583954	0.234584346	0.237980487	0.494220721	0	0	38130.99253	170132.9586	11.41355455	8.70E-05];

% paramMatrix <- each row := best fit params; ((row #) mod 3) := new seasonality function
paramMatrix=[6.046498	0.772609	0.465454	14.94935	0.984259	22.70853	18.91604	0	0	3.083237	1.06E-05; 11.81541	0.940728	0.50316	13.12023	0.971971	24.05734	10.78628	23843.62	2.529836	0.407502	9.83E-05; 8.970806	0.827267	0.48197	12.74008	0.002781	22.91748	19.21645	18.61354	54019.9	1.442842	4.85E-05;...
    6.629742318	0.497016921	0.758448622	13.888176	0.98017139	22.09773778	0.246581731	0	0	2.62136194	1.26E-06; 5.711284047	0.886456558	3.98E-01	11.92945054	9.81E-01	1.88E+01	6.75E-17	68905.59496	59119.62092	1.12E+01	6.09E-08; 14.91650856	0.083417469	0.972704762	14.97691685	0.000393962	23.6683787	0.205962629	70735.17542	65186.36826	2.449598712	8.36E-05; ...
    6.843815	0.135633	0.536052	10.18376	0.981556	2.675995	2.197288	0	0	10.85841	6.11E-05; 9.732335507	0.88375731	0.974724183	14.64355308	0.22065105	16.55353215	12.06242071	21603.46392	819.7085765	0.439535124	2.43E-05; 14.88164	0.633079	0.7087	14.11783	0.711675	2.650165	2.11021	14073.53	1.430811	0.7325	2.27E-05;...
    4.34273	0.727696	0.621459	1.288409	0.111454	0	0	0	0	10.28612	8.61E-05; 9.925002222	0.252880077	0.641201541	6.822712794	0.504776855	0	0	39253.58221	51046.7589	12.21180941	9.80E-08; 14.63670049	0.974515145	0.262123165	8.629760033	0.735186818	0	0	31819.01589	34597.09844	10.78641526	5.90E-05; ...
    14.01871	0.021676	0.437779	12.15198	0.978639	0.00E+00	0.00E+00	0	0	8.655468	2.43E-03];

% row order: (1) exp biting, (2) exp birth, (3) exp both, (4) cos bite,
% (5) cos birth, (6) cos both, (7) beta bite, (8) beta birth, (9) beta both,
% (10) lin precip bite, (11) lin precip birth, (12) lin precip both,
% (13) no seasonality


%Initial Values
IH0=infectiondata(1);       %initial infected humans
SH0=3.941e6-IH0;            %initial susceptible humans
IF0=0;                      %initial infected forest midges
SC0=(7*510/3)*(1/((426299/122.57)*(1.21+3.92)))*2054713;
SF0=SC0*x(10);              %initial susceptible forest midge <- PLACEHOLDER
IC0=0;                      %initial infected city midges
SA0=(1559256.365-11401)*1.7;%initial susceptible sloths
IA0=SA0*x(11);              %initial infected sloths  <- PLACEHOLDER
IHC0=cumulativedata(1);

initialvalues = [SA0 IA0 SF0 IF0 SH0 IH0 SC0 IC0 IHC0]'; % intial value matrix

fineMesh = linspace(1,tspan(end),1001); % for pretty pictures

%% looping through for subplot

nIDs = 12;
alphabet = ('a':'l').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')');
charForTitle = cellstr(charlbl);

figure(2)
for i=1:(height(paramMatrix)-1)
    x=paramMatrix(i,:);
    initialvalues(2) = SA0*x(11);
    initialvalues(3) = SC0*x(10);

    if mod(i,3) == 1
        forceType = [1 0 0];
    elseif mod(i,3) == 2
        forceType = [0 1 0];
    else
        forceType = [1 1 0];
    end

    if i <= 3
        forceType(3) = 1;
    elseif i<=6
        forceType(3) = 2;
    elseif i<=9
        forceType(3) = 3;
    else
        forceType(3) = 4;
    end

    %solve ode
    [~,y] = ode45(@(t,y) OROV_Model_auto(t,y,x,forceType),fineMesh,initialvalues);
    IH = y(:,6);
    
    [~,y2] = ode45(@(t,y) OROV_Model_auto(t,y,x,forceType),tspan,initialvalues);
    IH2 = y2(:,6);
    diff = IH2 - infectiondata';
    ep = sqrt(sum(diff.^2)/sum(infectiondata.^2)) * 100;


    subplot(4,3,i)
    p1 = plot(fineMesh,IH,'Color',colors(1,:),'LineWidth', 1.5);
    hold on
    p2 = scatter(tspan,infectiondata,24,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
    hold on
    p3 = scatter(tspan,infectiondata,24,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
    % legend('Simulated Data', 'Observed Data', 'Location', 'northwest')
    if i>9
        xlabel('Time (weeks)')
    end
    if (i==1) || (i==4) || (i==7) || (i==10)
        ylabel('Human Infections')
    end
    xlim([0 52])
    ylim([-5 1000])
    % legendEP = "EP = " + string(ep);
    legend([p1 p2 p3],{'Simulation', 'Observed', "EP = " + string(ep)},'Location','northwest')
    set(gca,'fontsize',12)
    title(charForTitle(i),'FontSize',14)
    % if i<12
    %     text(0.025,0.925,charlbl{i},'Units','normalized','FontSize',14)
    % else
    %     text(0.95,0.925,charlbl{i},'Units','normalized','FontSize',14)
    % end
    hold off
end

%% alternative to alternative to fig 3

nIDs = 5;
alphabet = ('a':'e').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')');
charForTitle = cellstr(charlbl);
pos = [2 3 5 6];
count=1;

figure(3)
for i=[4 8 1 11 13]
    x=paramMatrix(i,:);
    initialvalues(2) = SA0*x(11);
    initialvalues(3) = SC0*x(10);
    
    if mod(i,3) == 1
        forceType = [1 0 0];
    elseif mod(i,3) == 2
        forceType = [0 1 0];
    else
        forceType = [1 1 0];
    end

    if i <= 3
        forceType(3) = 1;
    elseif i<=6
        forceType(3) = 2;
    elseif i<=9
        forceType(3) = 3;
    else
        forceType(3) = 4;
    end

    if i==13
        forceType = [0 0 0];
    end

    %solve ode
    [~,y] = ode45(@(t,y) OROV_Model_auto(t,y,x,forceType),fineMesh,initialvalues);
    IH = y(:,6);
    
    [~,y2] = ode45(@(t,y) OROV_Model_auto(t,y,x,forceType),tspan,initialvalues);
    IH2 = y2(:,6);
    diff = IH2 - infectiondata';
    ep = sqrt(sum(diff.^2)/sum(infectiondata.^2)) * 100;

    if i~=4
        subplot(2,3,pos(count-1))
        p1 = plot(fineMesh,IH,'Color',colors(1,:),'LineWidth', 1.5);
        hold on
        p2 = scatter(tspan,infectiondata,24,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
        hold on
        p3 = plot(0, 0, "o", 'color', 'none', 'MarkerSize', 10);
        % p3 = stem(nan, nan,'filled', 'color', colors(1,:),'LineWidth',1);
    else
        subplot(2,3,[1 4])
        p1 = plot(fineMesh,IH,'Color',colors(1,:),'LineWidth', 1.5);
        hold on
        p2 = scatter(tspan,infectiondata,24,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
        hold on
        p3 = plot(0, 0, "o", 'color', 'none', 'MarkerSize', 10);
        % p3 = stem(nan, nan,'filled', 'color', colors(1,:),'LineWidth',1);
        xlabel('Time (weeks)')
        ylabel('Human Infections')
    end

    if count>1
        if pos(count-1) >= 4
        xlabel('Time (weeks)')
        end
        if (pos(count-1) == 2) || (pos(count-1) == 5)
            ylabel('Human Infections')
        end
    end
    xlim([0 52])
    ylim([-5 1000])
    % legendEP = "EP = " + string(ep);
    legend([p1 p2 p3],{'Simulation', 'Observed', "EP = " + string(ep)},'Location','northwest')
    % legend('Simulation','Observed', [newline "EP = " + string(ep)],'Location','northwest')
    set(gca,'fontsize',18)
    title(charForTitle(count),'FontSize',24)
    hold off
    count=count+1;
end