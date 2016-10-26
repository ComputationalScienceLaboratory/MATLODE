clear global;

% User Supplied Functions
Ode_Function        = @swe_Function_mex;

disp( 'Shallow Water Equations Problem: ' );
Tspan = [0 3.2];
y0 = swe_Initialize_y0();
NVAR = max(size(y0));

% SDIRK Coefficient     ERK Coefficient     RK Coefficient      Ros Coefficient
%   1: Sdirk2A              1: Erk23            1: Rada2A           1: Ros2
%   2: Sdirk2B              2: Erk3_Heun        2: Lobatto3C        2: Ros3
%   3: Sdirk3A              3: Erk43            3: Gauss            3: Ros4
%   4: Sdirk4A              4: Dopri5           4: Radau1A          4: Rodas3
%   5: Sdirk4B              5: Verme            5: Lobatto3A        5: Rodas4
%                           6: Dopri853

RelTol = ones(NVAR,1).*1e-10;
AbsTol = ones(NVAR,1).*1e-10;

yDimension = length(y0);

Options = MATLODE_OPTIONS( 'AbsTol',          AbsTol, ...
    'RelTol',          RelTol, ...
    'AbsTol_TLM',      ones(yDimension,yDimension).*AbsTol(1).*10, ...
    'RelTol_TLM',      ones(yDimension,yDimension).*RelTol(1).*10, ...
    'Y_TLM',           eye(yDimension,yDimension), ...
    'NTLM',            yDimension, ...
    'AbsTol_ADJ',      AbsTol.*10, ...
    'RelTol_ADJ',      RelTol.*10, ...
    'NADJ',            yDimension, ...
    'Desired_Mode',    -1*adj_mode, ...
    'storeCheckpoint', true, ...
    'displayStats',    false, ...
    'displaySteps',    true, ...
    'Hmin',            0, ...
    'Hmax',            0, ...
    'Hstart',          0, ...
    'FacMin',          0, ...
    'FacMax',          0, ...
    'FacRej',          0, ...
    'FacSafe',         0, ...
    'Qmin',            0, ...
    'Qmax',            0, ...
    'WarningConfig',   0, ...
    'Autonomous',      0, ... % 1
    'ITOL',            0, ... % 2
    'Method',          5, ... % 3
    'Max_no_steps',    0, ... % 4
    'NewtonMaxit',     0, ... % 5
    'StartNewton',     0, ... % 6
    'SaveLU',          0, ... % 8
    'TLMNewtonEst',    0, ... % 9
    'TLMtruncErr',     0, ... % 12
    'FDAprox',         0, ... % 13
    'AdjointSolve',    0, ... % 14
    'DirectADJ',       0, ... % 16
    'ChunkSize',       50 );

% Compute Reference
disp('Reference');
options = odeset('AbsTol', 1d-13,'RelTol',1d-13);
[ T_Ref, Y_Ref ] = ode15s(Ode_Function,Tspan,y0, options);

Npoints = 11;
TOLS = logspace(-7,-2,Npoints);

ISTATUS_Array = cell(Npoints,1);
RSTATUS_Array = cell(Npoints,1);

% Preallocate for speed
elapsedTime_ode45 = zeros(1,Npoints);
elapsedTime_ode113 = zeros(1,Npoints);
elapsedTime_ERK_FWD_1 = zeros(1,Npoints);
elapsedTime_ERK_FWD_4 = zeros(1,Npoints);
elapsedTime_ERK_FWD_5 = zeros(1,Npoints);

nsteps_ode45 = zeros(1,Npoints);
nsteps_ode113 = zeros(1,Npoints);
nsteps_ERK_FWD_1 = zeros(1,Npoints);
nsteps_ERK_FWD_4 = zeros(1,Npoints);
nsteps_ERK_FWD_5 = zeros(1,Npoints);

error_ode45 = zeros(1,Npoints);
error_ode113 = zeros(1,Npoints);
error_ERK_FWD_1 = zeros(1,Npoints);
error_ERK_FWD_4 = zeros(1,Npoints);
error_ERK_FWD_5 = zeros(1,Npoints);

for ipt=1:Npoints
    RelTol=TOLS(ipt);
    AbsTol = RelTol;
    Options = MATLODE_OPTIONS( Options, ...
        'displaySteps', false, ...
        'displayStats', false, ...
        'storeCheckpoint', false, ...
        'RelTol',     ones(yDimension,1)*RelTol, ...
        'AbsTol',     ones(yDimension,1)*AbsTol, ...
        'Max_no_steps', 400000, ...
        'WarningConfig', 0);
    Options_MATLAB = odeset('AbsTol',AbsTol,'RelTol',RelTol);
    
    disp('ode45');
    tic;
    output_ode45 = ode45(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode45(ipt) = toc;
    T = transpose(output_ode45.x);
    Y = transpose(output_ode45.y);
    nsteps_ode45(ipt) = output_ode45.stats.nsteps;
    error_ode45(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ode113');
    tic;
    output_ode113 = ode113(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode113(ipt) = toc;
    Y = transpose(output_ode113.y);
    nsteps_ode113(ipt) = output_ode113.stats.nsteps;
    error_ode113(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ERK_FWD (Erk23)');
    Options = MATLODE_OPTIONS(Options,'Method',1);
    [ ~, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_ERK_FWD_1(ipt) = Stats.RSTATUS.Etime;
    nsteps_ERK_FWD_1(ipt) = Stats.ISTATUS.Nstp;
    error_ERK_FWD_1(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ERK_FWD (Dorpi5)');
    Options = MATLODE_OPTIONS(Options,'Method',4);
    [ ~, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_ERK_FWD_4(ipt) = Stats.RSTATUS.Etime;
    nsteps_ERK_FWD_4(ipt) = Stats.ISTATUS.Nstp;
    error_ERK_FWD_4(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ERK_FWD (Dopri853)');
    Options = MATLODE_OPTIONS(Options,'Method',5);
    [ ~, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_ERK_FWD_5(ipt) = Stats.RSTATUS.Etime;
    nsteps_ERK_FWD_5(ipt) = Stats.ISTATUS.Nstp;
    error_ERK_FWD_5(ipt) = relativeError(Y(end,:),Y_Ref(end,:));

end

figure;
loglog(elapsedTime_ode45,error_ode45,'--o', ...
    elapsedTime_ode113,error_ode113, '--^', ...
    elapsedTime_ERK_FWD_1,error_ERK_FWD_1, '-x', ...
    elapsedTime_ERK_FWD_4,error_ERK_FWD_4, '-s', ...
    elapsedTime_ERK_FWD_5,error_ERK_FWD_5, '-d');
legend('ode45', ...
    'ode113', ...
    'MATLODE\_ERK\_FWD (Erk23)', ...
    'MATLODE\_ERK\_FWD (Dopri5)',...
    'MATLODE\_ERK\_FWD (Dopri853)');
xlabel('CPU Time (sec)');
ylabel('Relative Error');
title('CPU vs. Relative Error (Shallow Water Equations)');

figure;
loglog(nsteps_ode45,error_ode45,'--o', ...
    nsteps_ode113,error_ode113, '--^', ...
    nsteps_ERK_FWD_1,error_ERK_FWD_1, '-x', ...
    nsteps_ERK_FWD_4,error_ERK_FWD_4, '-s', ...
    nsteps_ERK_FWD_5,error_ERK_FWD_5, '-d');
legend('ode45', ...
    'ode113', ...
    'MATLODE\_ERK\_FWD (Erk23)', ...
    'MATLODE\_ERK\_FWD (Dopri5)',...
    'MATLODE\_ERK\_FWD (Dopri853)');
xlabel('Number of Steps');
ylabel('Relative Error');
title('Number of Steps vs. Relative Error (Shallow Water Equations)');