resultsPath = './Lorenz_Results';

Ode_Function        = @lorenz_Function;
Ode_Jacobian        = @lorenz_Jacobian;

Tspan = [0 20];
y0 = [5;5;17];

RelTol = ones(3,1)*1e-8;
AbsTol = ones(3,1)*1e-8;
         
RelTol_TLM = ones(3,1)*1e-7;
AbsTol_TLM = ones(3,1)*1e-7;

Y_TLM = zeros(3,3);
for k=1:3
    Y_TLM(k,k) = 1;
end

Options = MATLODE_OPTIONS( 'AbsTol',          AbsTol, ...
                           'RelTol',          RelTol, ...
                           'Jacobian',        Ode_Jacobian, ...
                           'AbsTol_TLM',      AbsTol_TLM, ...
                           'RelTol_TLM',      RelTol_TLM, ...
                           'Y_TLM',           Y_TLM, ...
                           'storeCheckpoint', true, ...
                           'NTLM',            3);
                                              
disp( 'Solving problem with SDIRK_TLM_Integrator: ');
[ T, Y, Sens, Stats ] = MATLODE_SDIRK_TLM_Integrator( Ode_Function, Tspan, y0, Options );

labels = { 'beta', 'rho', 'sigma' };

plot3( Y(:,1), Y(:,2), Y(:,3) );
title('Lorenz', 'FontSize', 20);
xlabel('beta', 'FontSize', 20);
ylabel('rho', 'FontSize', 20);
zlabel('sigma', 'FontSize', 20);
saveas(gcf,strcat(resultsPath,'/LORENZ_Solution.pdf'),'pdf');
saveas(gcf,strcat(resultsPath,'/LORENZ_Solution.fig'),'fig');

% BETA
[~, new_indices] = sort(abs(Sens(:,1))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(new_indices,1));
title('Y1');
set(gca,'YTickLabel',sorted_labels);
saveas(gcf,strcat(resultsPath,'/LORENZ_BETA_SDIRK_TLM.pdf'),'pdf');
saveas(gcf,strcat(resultsPath,'/LORENZ_BETA_SDIRK_TLM.fig'),'fig');

% RHO
[~, new_indices] = sort(abs(Sens(:,2))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(new_indices,2));
title('Y2');
set(gca,'YTickLabel',sorted_labels);
saveas(gcf,strcat(resultsPath,'/LORENZ_RHO_SDIRK_TLM.pdf'),'pdf');
saveas(gcf,strcat(resultsPath,'/LORENZ_RHO_SDIRK_TLM.fig'),'fig');

% SIGMA 
[~, new_indices] = sort(abs(Sens(:,3))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(new_indices,3));
title('Y3');
set(gca,'YTickLabel',sorted_labels);
saveas(gcf,strcat(resultsPath,'/LORENZ_SIGMA_SDIRK_TLM.pdf'),'pdf');
saveas(gcf,strcat(resultsPath,'/LORENZ_SIGMA_SDIRK_TLM.fig'),'fig');