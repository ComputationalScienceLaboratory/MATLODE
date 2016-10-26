resultsPath = './CBM4_Results';

cbm_Parameters;
cbm_Global_defs;
cbm_Sparse;
cbm_Monitor;
cbm_JacobianSP;
cbm_HessianSP;
cbm_StoichiomSP;

TSTART = 12*3600;
TEND = TSTART + 7*24*3600;

DT = 60.;
TEMP = 298;

cbm_Initialize;

Ode_Function        = @cbm_Fun_Chem;
Ode_Jacobian        = @cbm_Jac_Chem;

TIME = TSTART;

Tspan = [TSTART TEND];

RelTol = ones(32,1)*1e-8;
AbsTol = ones(32,1)*1e-8;
         
RelTol_ADJ = ones(32,1)*1e-7;
AbsTol_ADJ = ones(32,1)*1e-7;

y0 = VAR;

Lambda = zeros(32,81);
for k=1:81
    if ( k > 32 )
        break;
    end
    Lambda(k,k) = 1;
end

Options = MATLODE_OPTIONS( 'AbsTol',          AbsTol, ...
                           'RelTol',          RelTol, ...
                           'Jacobian',        Ode_Jacobian, ...
                           'AbsTol_ADJ',      AbsTol_ADJ, ...
                           'RelTol_ADJ',      RelTol_ADJ, ...
                           'NADJ',            81, ...
                           'storeCheckpoint', true, ...
                           'displayStats',    true, ...
                           'displaySteps',    false, ...
                           'Lambda',          Lambda, ...
                           'Max_no_steps',    400000, ...
                           'ChunkSize',       50);

                       
disp( 'Solving problem with SDIRK_ADJ_Integrator: ');
[ T, Y, Sens, Stats ] = MATLODE_SDIRK_ADJ_Integrator( Ode_Function, Tspan, y0, Options );

labels = { '1', 'O+O2+M=O3', '3', '4', '5', '6', '7', '8', '9', '10', ...
    '11', 'O3+OH=HO2', '13', '14', '15', 'NO3+NO2=NO+NO2', '17', '18', '19', '20', ...
    'NO+NO2+H2O=2HONO', '22', '23', '24', ' 2HONO=NO+NO2', 'OH+NO2=HNO3', '27', '28', '29', '30', ...
    'OH+PNA=NO2', '32', '33', '34', '35', '36', '37', '38', '39', '40', ...
    '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', ...
    '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', ...
    '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', ...
    '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', ...
    '81'};

% N2O5
[~, new_indices] = sort(abs(Sens(6,:))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(6,new_indices(75:end)));
title('N2O5');
set(gca,'YTickLabel',sorted_labels(75:end));
saveas(gcf,strcat(resultsPath,'/CBM4_N2O5_SDIRK_Adjoint.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_N2O5_SDIRK_Adjoint.fig'),'fig');

% HONO
[~, new_indices] = sort(abs(Sens(9,:))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(9,new_indices(75:end)));
title('HONO');
set(gca,'YTickLabel',sorted_labels(75:end));
saveas(gcf,strcat(resultsPath,'/CBM4_HONO_SDIRK_Adjoint.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_HONO_SDIRK_Adjoint.fig'),'fig');

% HNO3 
[~, new_indices] = sort(abs(Sens(12,:))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(12,new_indices(75:end)));
title('HNO3');
set(gca,'YTickLabel',sorted_labels(75:end));
saveas(gcf,strcat(resultsPath,'/CBM4_HNO3_SDIRK_Adjoint.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_HNO3_SDIRK_Adjoint.fig'),'fig');

% O3
[~, new_indices] = sort(abs(Sens(25,:))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(25,new_indices(75:end)));
title('O3');
set(gca,'YTickLabel',sorted_labels(75:end));
saveas(gcf,strcat(resultsPath,'/CBM4_O3_SDIRK_Adjoint.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_O3_SDIRK_Adjoint.fig'),'fig');

% NO2 
[~, new_indices] = sort(abs(Sens(26,:))); % sorts in *ascending* order
sorted_labels = labels(new_indices);
figure();
barh(Sens(26,new_indices(75:end)));
title('NO2');
set(gca,'YTickLabel',sorted_labels(75:end));
saveas(gcf,strcat(resultsPath,'/CBM4_NO2_SDIRK_Adjoint.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_NO2_SDIRK_Adjoint.fig'),'fig');