clear global;

% User Supplied Functions
Ode_Function        = @cbm_Fun_Chem;
Ode_Jacobian        = @cbm_Jac_Chem;

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

C(1:32) = VAR(1:32);
C((32+1):33) = FIX(1:1);

TIME = TSTART;
Tspan = [TSTART TEND];

RelTol = ones(32,1)*1e-3;
AbsTol = ones(32,1)*1e+3;

y0 = VAR;
yDimension = length(VAR);

Options = MATLODE_OPTIONS( 'AbsTol',          AbsTol, ...
    'RelTol',          RelTol, ...
    'Jacobian',        Ode_Jacobian, ...
    'AbsTol_TLM',      AbsTol.*10, ...
    'RelTol_TLM',      RelTol.*10, ...
    'Y_TLM',           eye(yDimension,yDimension), ...
    'NTLM',            yDimension, ...
    'AbsTol_ADJ',      AbsTol.*10, ...
    'RelTol_ADJ',      RelTol.*10, ...
    'NADJ',            yDimension, ...
    'storeCheckpoint', false, ...
    'displayStats',    true, ...
    'displaySteps',    false, ...
    'ChunkSize',       50 );


% Compute Reference
disp('Reference ode15s');
options = odeset('AbsTol', 1d-8,'RelTol',1d-8,'Jacobian',Ode_Jacobian);
[ ~, Y_Ref ] = ode15s(Ode_Function,Tspan,y0, options);

disp('Reference ode23s');
options = odeset('AbsTol', 1d-8,'RelTol',1d-8,'Jacobian',Ode_Jacobian);
[ ~, Y_Ref_ode23s ] = ode23s(Ode_Function,Tspan,y0, options);

Npoints = 6;
TOLS = logspace(-6,-1,Npoints);

ISTATUS_Array = cell(Npoints,1);
RSTATUS_Array = cell(Npoints,1);

% Used for constructing legend and method
rosMethodNumber = 5;
rkMethodNumber = 4;
sdirkMethodNumber = 5;

switch (rosMethodNumber)
    case 1
        rosMethodName = 'Ros2';
    case 2
        rosMethodName = 'Ros3';
    case 3
        rosMethodName = 'Ros4';
    case 4
        rosMethodName = 'Rodas3';
    case 5
        rosMethodName = 'Rodas4';
    otherwise
        rosMethodName = 'Ros4';
end

switch (sdirkMethodNumber)
    case 1
        sdirkMethodName = 'Sdirk2A';
    case 2
        sdirkMethodName = 'Sdirk2B';
    case 3
        sdirkMethodName = 'Sdirk3A';
    case 4
        sdirkMethodName = 'Sdirk3B';
    case 5
        sdirkMethodName = 'Sdirk4A';
    case 6
        sdirkMethodName = 'Sdirk4B';
    otherwise
        sdiirkMethodName = 'Sdirk4A';
end

switch (rkMethodNumber)
    case 1
        rkMethodName = 'Radau2A';
    case 2
        rkMethodName = 'Lobatto3C';
    case 3
        rkMethodName = 'Gauss';
    case 4
        rkMethodName = 'Radau1A';
    case 5
        rkMethodName = 'Lobatto3A';
    otherwise
        rkMethodName = 'Lobatto3A';
end


% Preallocate for speed
elapsedTime_ode45 = zeros(1,Npoints);
elapsedTime_ode15s = zeros(1,Npoints);
elapsedTime_ode23 = zeros(1,Npoints);
elapsedTime_ode113 = zeros(1,Npoints);
elapsedTime_ode23s = zeros(1,Npoints);
elapsedTime_ode23t = zeros(1,Npoints);
elapsedTime_ode23tb = zeros(1,Npoints);
elapsedTime_ERK_FWD_1 = zeros(1,Npoints);
elapsedTime_ERK_FWD_4 = zeros(1,Npoints);
elapsedTime_ERK_FWD_5 = zeros(1,Npoints);
elapsedTime_ROS_FWD = zeros(1,Npoints);
elapsedTime_RK_FWD = zeros(1,Npoints);
elapsedTime_SDIRK_FWD = zeros(1,Npoints);

nsteps_ode45 = zeros(1,Npoints);
nsteps_ode15s = zeros(1,Npoints);
nsteps_ode23 = zeros(1,Npoints);
nsteps_ode113 = zeros(1,Npoints);
nsteps_ode23s = zeros(1,Npoints);
nsteps_ode23t = zeros(1,Npoints);
nsteps_ode23tb = zeros(1,Npoints);
nsteps_ERK_FWD_1 = zeros(1,Npoints);
nsteps_ERK_FWD_4 = zeros(1,Npoints);
nsteps_ERK_FWD_5 = zeros(1,Npoints);
nsteps_ROS_FWD = zeros(1,Npoints);
nsteps_RK_FWD = zeros(1,Npoints);
nsteps_SDIRK_FWD = zeros(1,Npoints);

error_ode45 = zeros(1,Npoints);
error_ode15s = zeros(1,Npoints);
error_ode23 = zeros(1,Npoints);
error_ode113 = zeros(1,Npoints);
error_ode23s = zeros(1,Npoints);
error_ode23t = zeros(1,Npoints);
error_ode23tb = zeros(1,Npoints);
error_ERK_FWD_1 = zeros(1,Npoints);
error_ERK_FWD_4 = zeros(1,Npoints);
error_ERK_FWD_5 = zeros(1,Npoints);
error_ROS_FWD = zeros(1,Npoints);
error_RK_FWD = zeros(1,Npoints);
error_SDIRK_FWD = zeros(1,Npoints);

for ipt=1:Npoints
    RelTol=TOLS(ipt);
    AbsTol = RelTol;
    Options = MATLODE_OPTIONS( Options, ...
        'displaySteps', false, ...
        'displayStats', false, ...
        'storeCheckpoint', false, ...
        'RelTol',     ones(yDimension,1)*RelTol, ...
        'AbsTol',     ones(yDimension,1)*AbsTol, ...
        'Jacobian',   Ode_Jacobian, ...
        'Max_no_steps', 400000, ...
        'WarningConfig', 0);
    Options_MATLAB = odeset('AbsTol',AbsTol,'RelTol',RelTol,'Jacobian',Ode_Jacobian);
    
    disp('ode15s');
    tic;
    output_ode15s = ode15s(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode15s(ipt) = toc;
    T = transpose(output_ode15s.x);
    Y = transpose(output_ode15s.y);
    nsteps_ode15s(ipt) = output_ode15s.stats.nsteps;
    error_ode15s(ipt) = relativeError(Y(end,:),Y_Ref_ode23s(end,:));
    
    disp('ode23s');
    tic;
    output_ode23s = ode23s(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode23s(ipt) = toc;
    T = transpose(output_ode23s.x);
    Y = transpose(output_ode23s.y);
    nsteps_ode23s(ipt) = output_ode23s.stats.nsteps;
    error_ode23s(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ode23t');
    tic;
    output_ode23t = ode23t(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode23t(ipt) = toc;
    T = transpose(output_ode23t.x);
    Y = transpose(output_ode23t.y);
    nsteps_ode23t(ipt) = output_ode23t.stats.nsteps;
    error_ode23t(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ode23tb');
    tic;
    output_ode23tb = ode23tb(Ode_Function,Tspan,y0,Options_MATLAB);
    elapsedTime_ode23tb(ipt) = toc;
    T = transpose(output_ode23tb.x);
    Y = transpose(output_ode23tb.y);
    nsteps_ode23tb(ipt) = output_ode23tb.stats.nsteps;
    error_ode23tb(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('ROS_FWD');
    % Method: Ros2 {1}   | Ros3 {2}   | Ros4 {3} |
    %         Rodas3 {4} | Rodas4 {5}
    Options = MATLODE_OPTIONS(Options, 'Method', rosMethodName);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ T, Y, Stats ] = MATLODE_ROS_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_ROS_FWD(ipt) = Stats.RSTATUS.Etime;
    nsteps_ROS_FWD(ipt) = Stats.ISTATUS.Nstp;
    error_ROS_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('RK_FWD');
    Options = MATLODE_OPTIONS(Options, 'Method', rkMethodName);
    [ T, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('SDIRK_FWD');
    Options = MATLODE_OPTIONS(Options, 'Method', sdirkMethodName);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ T, Y, Stats ] = MATLODE_SDIRK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_SDIRK_FWD(ipt) = Stats.RSTATUS.Etime;
    nsteps_SDIRK_FWD(ipt) = Stats.ISTATUS.Nstp;
    error_SDIRK_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
end

rosTitle = ['MATLODE\_ROS\_FWD (', rosMethodName, ')'];
sdirkTitle = ['MATLODE\_SDIRK\_FWD (', sdirkMethodName, ')'];
rkTitle = ['MATLODE\_RK\_FWD (', rkMethodName, ')'];

figure;
loglog(elapsedTime_ode15s,error_ode15s,'--+', ...
    elapsedTime_ode23s,error_ode23s, '--x', ...
    elapsedTime_ode23t,error_ode23t, '--s', ...
    elapsedTime_ode23tb,error_ode23tb, '--d', ...
    elapsedTime_ROS_FWD,error_ROS_FWD, '-+', ...
    elapsedTime_RK_FWD,error_RK_FWD, '-*', ...
    elapsedTime_SDIRK_FWD,error_SDIRK_FWD, '-^');
CPU_Legend = legend('ode15s', ...
    'ode23s', ...
    'ode23t', ...
    'ode23tb', ...
    rosTitle, ...
    rkTitle, ...
    sdirkTitle);
xlabel('CPU Time (sec)');
ylabel('Relative Error');
CPU_Title = title('CPU vs. Relative Error (CBM4)');
set(gca,'FontSize',16);
set(CPU_Title,'FontSize',16);
set(CPU_Legend,'FontSize',16);
axis([5e-1 10e2 10e-10 10e-2]);
saveas(gcf,'CBM4_Results/CBM4_Implicit_CpuResults.eps','epsc');
saveas(gcf,'CBM4_Results/CBM4_Implicit_CpuResults.fig','fig');
saveas(gcf,'CBM4_Results/CBM4_Implicit_CpuResults.pdf','pdf');

figure;
loglog(nsteps_ode15s,error_ode15s,'--+', ...
    nsteps_ode23s,error_ode23s, '--x', ...
    nsteps_ode23t,error_ode23t, '--s', ...
    nsteps_ode23tb,error_ode23tb, '--d', ...
    nsteps_ROS_FWD,error_ROS_FWD, '-+', ...
    nsteps_RK_FWD,error_RK_FWD, '-*', ...
    nsteps_SDIRK_FWD,error_SDIRK_FWD, '-^');
NSteps_Legend = legend('ode15s', ...
    'ode23s', ...
    'ode23t', ...
    'ode23tb', ...
    rosTitle', ...
    rkTitle, ...
    sdirkTitle);
xlabel('Number of Steps');
ylabel('Relative Error');
NSteps_Title = title('Number of Steps vs. Relative Error (CBM4)');
set(gca,'FontSize',16);
set(NSteps_Title,'FontSize',16);
set(NSteps_Legend,'FontSize',16);
axis([10e1 10e6 10e-10 10e-1]);
saveas(gcf,'CBM4_Results/CBM4_Implicit_NstepsResults.eps','epsc');
saveas(gcf,'CBM4_Results/CBM4_Implicit_NstepsResults.fig','fig');
saveas(gcf,'CBM4_Results/CBM4_Implicit_NstepsResults.pdf','pdf');

