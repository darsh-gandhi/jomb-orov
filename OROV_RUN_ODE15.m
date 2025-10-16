%This is the function which sets up the ODE problem and also what we are minimizing

function value=OROV_RUN_ODE15(z,data,new_initialvalues,tspan,w,nu)
        [t,y] = ode15s(@(t,y) OROV_Model(t,y,z),tspan,new_initialvalues);
        SA = y(:,1);
        IA = y(:,2);
        SF = y(:,3);
        IF = y(:,4);
        SH = y(:,5);
        IH = y(:,6);
        SC = y(:,7);
        IC = y(:,8);
        IHC = y(:,9);
       
%Difference Data and Model
% diff = expected infection data vs actual infection data

% Weekly data
diff = IH - data';

% % Cumulative data
% diff = IHC - data';

%Cost functional
value = w*0.5*sum(sum((diff).^2)) + nu*0.5*sum(z.^2);
end