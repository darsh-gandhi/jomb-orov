close all; clear all;

paramNames = {'$b$', '$\pi_{VH}$', '$\pi_{HV}$', '$k$', '$p$', '$t^*$', '$a$', '$\Lambda_F$', '$\Lambda_C$', '$r_F$', '$r_C$'};
baseError = 11.94;
errors = [14.65; 12.77; 12.59; 11.98; 35.91; 13.04; 12.38; 12.49; 11.94; 11.99; 11.96];

diff = abs(errors - baseError);

[errorSorted, idx] = sort(diff, 'descend');
paramNamesSorted = paramNames(idx);

figure;
barh(errorSorted, 'FaceColor',[0.2 0.6 0.8]);
set(gca, 'ytick', 1:length(paramNamesSorted),'yticklabel', paramNamesSorted,'TickLabelInterpreter', 'latex','YDir', 'reverse','FontSize',14); 
xlabel('Change in EP (%)','FontSize',14);
set(gca, 'FontSize', 14);
