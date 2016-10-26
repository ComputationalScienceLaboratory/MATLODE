clear global;

% error: -1: relative root mean squared (rRMS)
%        -2: relative error
%        -3: absolute error
error_mode = -1;

% integrator: -1: ERK_FWD_Integrator
%             -2: ROS_FWD_Integrator
%             -3: RK_FWD_Integrator
%             -4: SDIRK_FWD_Integrator
%
%             -5: ERK_TLM_Integrator
%             -6: ROS_TLM_Integrator
%             -7: RK_TLM_Integrator
%             -8: SDIRK_TLM_Integrator
%
%              -9: ERK_ADJ_Integrator
%             -10: ROS_ADJ_Integrator
%             -11: RK_ADJ_Integrator
%             -12: SDIRK_ADJ_Integrator              
integrator = -4;

% Adjoint: -1: xxxxADJ1 w/ Quadrature
%          -2: xxxxADJ1 no Quadrature
%          -3: xxxxADJ2 w/ Quadrature
%          -4: xxxxADJ2 no Quadrature
adj_mode = -4;

% User Supplied Functions
Ode_Function        = @cbm_Fun_Chem;
Ode_Jacobian        = @cbm_Jac_Chem;
% Ode_Hess_vec        = @cbm4_Hess_vec_mex;
% Ode_Hesstr_vec      = @cbm4_Hesstr_vec_mex;
% Ode_Jacp            = @cbm4_Jacp;
% Ode_Hesstr_vec_f_py = @cbm4_Hesstr_vec_f_py;
% Ode_QFun            = @cbm4_QFun;
% Ode_Mu              = @cbm4_Mu;
% Ode_Quadrature      = @cbm4_Quadrature;
% Ode_DRDP            = @cbm4_DRDP;
% Ode_DRDY            = @cbm4_DRDY;
% Ode_Hesstr_vec_r_py = @cbm4_Hesstr_vec_r_py;
% Ode_Hesstr_vec_r    = @cbm4_Hesstr_vec_r;

resultsPath = './CBM4_Results';

%profile clear;
%profile on;

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

% SDIRK Coefficient     ERK Coefficient     RK Coefficient      Ros Coefficient
%   1: Sdirk2A              1: Erk23            1: Rada2A           1: Ros2
%   2: Sdirk2B              2: Erk3_Heun        2: Lobatto3C        2: Ros3
%   3: Sdirk3A              3: Erk43            3: Gauss            3: Ros4
%   4: Sdirk4A              4: Dopri5           4: Radau1A          4: Rodas3
%   5: Sdirk4B              5: Verme            5: Lobatto3A        5: Rodas4
%                           6: Dopri853

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
                           'Desired_Mode',    -1*adj_mode, ...
                           'storeCheckpoint', true, ...
                           'displayStats',    true, ...
                           'displaySteps',    false, ...
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
                           'Method',          0, ... % 3
                           'Max_no_steps',    0, ... % 4
                           'NewtonMaxit',     0, ... % 5
                           'StartNewton',     0, ... % 6
                           'DirectTLM',       0, ... % 7  % FATODE Needs this update
                           'SaveLU',          0, ... % 8                        
                           'TLMNewtonEst',    0, ... % 9
                           'TLMtruncErr',     0, ... % 12
                           'FDAprox',         0, ... % 13                       
                           'AdjointSolve',    0, ... % 14
                           'DirectADJ',       0, ... % 16                           
                           'ChunkSize',       50 );
                                     
if ( integrator == -5 )
    Options = MATLODE_OPTIONS( Options, ...
                               'FDIncrement', 0, ...
                               'FDAprox',     0, ...
                               'NTLM',        yDimension );
end
                       
if ( integrator == -6 ) 
    Options = MATLODE_OPTIONS( Options, ...
                                  'Hess_vec', Ode_Hess_vec );
end
         

if ( integrator == -7 )
    Options = MATLODE_OPTIONS( Options, ...
                               'FDIncrement', 10^-6, ...
                               'FDAprox',     1, ...
                               'NTLM',        yDimension );
end

if ( integrator == -10 )
    switch ( adj_mode )
        case -1
            disp( 'RosenbrockADJ1 w/ Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'NP', 1, ...
                                          'Hesstr_vec', Ode_Hesstr_vec, ...
                                          'Jacp', Ode_Jacp, ...
                                          'Lambda', eye(yDimension), ...
                                          'Hesstr_vec_f_py', Ode_Hesstr_vec_f_py, ...
                                          'QFun', Ode_QFun, ...
                                          'Mu', Ode_Mu, ...
                                          'Quadrature', Ode_Quadrature, ...
                                          'DRDP', Ode_DRDP, ...
                                          'DRDY', Ode_DRDY, ...
                                          'Hesstr_vec_r_py', Ode_Hesstr_vec_r_py, ...
                                          'Hesstr_vec_r', Ode_Hesstr_vec_r );
        case -2
            disp( 'RosenbrockADJ1 no Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'Hesstr_vec', Ode_Hesstr_vec, ...
                                          'NP', 1, ...
                                          'Jacp', Ode_Jacp, ...
                                          'Lambda', eye(yDimension), ...
                                          'Hesstr_vec_f_py', Ode_Hesstr_vec_f_py, ...
                                          'Mu', Ode_Mu ); 
        case -3
            disp( 'RosenbrockADJ2 w/ Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'Hesstr_vec', Ode_Hesstr_vec, ...
                                          'Lambda', eye(yDimension), ...
                                          'QFun', Ode_QFun, ...
                                          'Quadrature', Ode_Quadrature, ...
                                          'DRDY', Ode_DRDY, ...
                                          'Hesstr_vec_r', Ode_Hesstr_vec_r );
        case -4
            disp( 'RosenbrockADJ2 no Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'Lambda', eye(yDimension), ...
                                          'Hesstr_vec', Ode_Hesstr_vec );
        otherwise
            disp('Error: Choose Rosenbrock Adjoint Mode');
    end
end

if ( (integrator == -9) || (integrator == -11) || (integrator == -12) )
    switch ( adj_mode )
        case -1
            disp( 'ERK/RK/SDIRK ADJ1 w/ Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'NP', 1, ...
                                          'Jacp', Ode_Jacp, ...
                                          'Lambda', eye(yDimension), ...
                                          'QFun', Ode_QFun, ...
                                          'Mu', Ode_Mu, ...
                                          'Quadrature', Ode_Quadrature, ...
                                          'DRDP', Ode_DRDP, ...
                                          'DRDY', Ode_DRDY );
        case -2
            disp( 'ERK/RK/SDIRK ADJ1 no Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'NP', 1, ...
                                          'Jacp', Ode_Jacp, ...
                                          'Lambda', eye(yDimension), ...
                                          'Mu', Ode_Mu );
        case -3
            disp( 'ERK/RK/SDIRK ADJ2 w/ Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                          'NP', 1, ...
                                          'Lambda', eye(yDimension), ...
                                          'QFun', Ode_QFun, ...
                                          'Quadrature', Ode_Quadrature, ...
                                          'DRDY', Ode_DRDY );
        case -4
            disp( 'ERK/RK/SDIRK ADJ2 no Quadrature' );
            Options = MATLODE_OPTIONS( Options, ...
                                        'NP', 1, ...
                                        'Lambda', eye(yDimension) );
        otherwise
    end
end


% Compute Reference
disp('Reference');
options = odeset('AbsTol', 1d-8,'RelTol',1d-8,'Jacobian',Ode_Jacobian);
[ T_Ref, Y_Ref ] = ode15s(Ode_Function,Tspan,y0, options);

Npoints = 6;
TOLS = logspace(-6,-1,Npoints);

ISTATUS_Array = cell(Npoints,1);
RSTATUS_Array = cell(Npoints,1);


% Preallocate for speed
elapsedTime_RK_FWD_Radau2A = zeros(1,Npoints);
elapsedTime_RK_FWD_Radau2A_LU = zeros(1,Npoints);
elapsedTime_RK_FWD_Lobatto3C = zeros(1,Npoints);
elapsedTime_RK_FWD_Lobatto3C_LU = zeros(1,Npoints);
elapsedTime_RK_FWD_Gauss = zeros(1,Npoints);
elapsedTime_RK_FWD_Gauss_LU = zeros(1,Npoints);
elapsedTime_RK_FWD_Radau1A = zeros(1,Npoints);
elapsedTime_RK_FWD_Radau1A_LU = zeros(1,Npoints);

nsteps_RK_FWD_Radau2A = zeros(1,Npoints);
nsteps_RK_FWD_Radau2A_LU = zeros(1,Npoints);
nsteps_RK_FWD_Lobatto3C = zeros(1,Npoints);
nsteps_RK_FWD_Lobatto3C_LU = zeros(1,Npoints);
nsteps_RK_FWD_Gauss = zeros(1,Npoints);
nsteps_RK_FWD_Gauss_LU = zeros(1,Npoints);
nsteps_RK_FWD_Radau1A = zeros(1,Npoints);
nsteps_RK_FWD_Radau1A_LU = zeros(1,Npoints);

error_RK_FWD_Radau2A = zeros(1,Npoints);
error_RK_FWD_Radau2A_LU = zeros(1,Npoints);
error_RK_FWD_Lobatto3C = zeros(1,Npoints);
error_RK_FWD_Lobatto3C_LU = zeros(1,Npoints);
error_RK_FWD_Gauss = zeros(1,Npoints);
error_RK_FWD_Gauss_LU = zeros(1,Npoints);
error_RK_FWD_Radau1A = zeros(1,Npoints);
error_RK_FWD_Radau1A_LU = zeros(1,Npoints);

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
    
    disp('RK_FWD_Radau2A');
    Options = MATLODE_OPTIONS(Options, 'Method', 1, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 0);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Radau2A(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Radau2A(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Radau2A(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
        
    disp('RK_FWD_Radau2A_LU');
    Options = MATLODE_OPTIONS(Options, 'Method', 1, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Radau2A_LU(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Radau2A_LU(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Radau2A_LU(ipt) = relativeError(Y(end,:),Y_Ref(end,:));    
    
    disp('RK_FWD_Lobatto3C');
    Options = MATLODE_OPTIONS(Options, 'Method', 2, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 0);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Lobatto3C(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Lobatto3C(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Lobatto3C(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
        
    disp('RK_FWD_Lobatto3C_LU');
    Options = MATLODE_OPTIONS(Options, 'Method', 2, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Lobatto3C_LU(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Lobatto3C_LU(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Lobatto3C_LU(ipt) = relativeError(Y(end,:),Y_Ref(end,:));    
    
    disp('RK_FWD_Gauss');
    Options = MATLODE_OPTIONS(Options, 'Method', 3, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 0);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Gauss(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Gauss(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Gauss(ipt) = relativeError(Y(end,:),Y_Ref(end,:));  
    
    disp('RK_FWD_Gauss_LU');
    Options = MATLODE_OPTIONS(Options, 'Method', 3, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Gauss_LU(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Gauss_LU(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Gauss_LU(ipt) = relativeError(Y(end,:),Y_Ref(end,:));      
    
    disp('RK_FWD_Radau1A');
    Options = MATLODE_OPTIONS(Options, 'Method', 4, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 0);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Radau1A(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Radau1A(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Radau1A(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
    
    disp('RK_FWD_Radau1A_LU');
    Options = MATLODE_OPTIONS(Options, 'Method', 4, 'SdirkError', 1);
    Options = MATLODE_OPTIONS(Options, 'LU', 1);
    [ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );
    elapsedTime_RK_FWD_Radau1A_LU(ipt) = Stats.RSTATUS.Etime;
    nsteps_RK_FWD_Radau1A_LU(ipt) = Stats.ISTATUS.Nstp;
    error_RK_FWD_Radau1A_LU(ipt) = relativeError(Y(end,:),Y_Ref(end,:));    
    
end

% SDIRK Coefficient     ERK Coefficient     RK Coefficient      Ros Coefficient
%   1: Sdirk2A              1: Erk23            1: Rada2A           1: Ros2
%   2: Sdirk2B              2: Erk3_Heun        2: Lobatto3C        2: Ros3
%   3: Sdirk3A              3: Erk43            3: Gauss            3: Ros4
%   4: Sdirk4A              4: Dopri5           4: Radau1A          4: Rodas3
%   5: Sdirk4B              5: Verme            5: Lobatto3A        5: Rodas4

figure;
loglog(elapsedTime_RK_FWD_Radau2A,error_RK_FWD_Radau2A,'-+', ... 
    elapsedTime_RK_FWD_Radau2A_LU,error_RK_FWD_Radau2A_LU,'--+', ...
    elapsedTime_RK_FWD_Lobatto3C,error_RK_FWD_Lobatto3C, '-x', ...
    elapsedTime_RK_FWD_Lobatto3C_LU,error_RK_FWD_Lobatto3C_LU, '--x', ...    
    elapsedTime_RK_FWD_Gauss,error_RK_FWD_Gauss, '-s', ...
    elapsedTime_RK_FWD_Gauss_LU,error_RK_FWD_Gauss_LU, '--s', ...
    elapsedTime_RK_FWD_Radau1A,error_RK_FWD_Radau1A, '-d', ...
    elapsedTime_RK_FWD_Radau1A_LU,error_RK_FWD_Radau1A_LU, '--d');
legend('MATLODE\_RK\_FWD (Radau2A)', ...
    'MATLODE\_RK\_FWD (Radau2A LU)', ...
    'MATLODE\_RK\_FWD (Lobatto3C)', ...
    'MATLODE\_RK\_FWD (Lobatto3C LU)', ...
    'MATLODE\_RK\_FWD (Gauss)', ...
    'MATLODE\_RK\_FWD (Gauss LU)', ...
    'MATLODE\_RK\_FWD (Radau1A)', ...
    'MATLODE\_RK\_FWD (Radau1A LU)');
xlabel('CPU Time (sec)');
ylabel('Relative Error');
title('CPU vs. Relative Error (CBM4)');
saveas(gcf,strcat(resultsPath,'/CBM4_RK_ALL_Implicit_CpuResults.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_RK_ALL_Implicit_CpuResults.fig'),'fig');

figure;
loglog(nsteps_RK_FWD_Radau2A,error_RK_FWD_Radau2A,'-+', ...
    nsteps_RK_FWD_Radau2A_LU,error_RK_FWD_Radau2A_LU,'--+', ...
    nsteps_RK_FWD_Lobatto3C,error_RK_FWD_Lobatto3C, '-x', ...
    nsteps_RK_FWD_Lobatto3C_LU,error_RK_FWD_Lobatto3C_LU, '--x', ...
    nsteps_RK_FWD_Gauss,error_RK_FWD_Gauss, '-s', ...
    nsteps_RK_FWD_Gauss_LU,error_RK_FWD_Gauss_LU, '--s', ...
    nsteps_RK_FWD_Radau1A,error_RK_FWD_Radau1A, '-d', ...
    nsteps_RK_FWD_Radau1A_LU,error_RK_FWD_Radau1A_LU', '--d');
legend('MATLODE\_RK\_FWD (Radau2A)', ...
    'MATLODE\_RK\_FWD (Radau2A LU)', ...
    'MATLODE\_RK\_FWD (Lobatto3C)', ...
    'MATLODE\_RK\_FWD (Lobatto3C LU)', ...
    'MATLODE\_RK\_FWD (Gauss)', ...
    'MATLODE\_RK\_FWD (Gauss LU)', ...
    'MATLODE\_RK\_FWD (Radau1A)', ...
    'MATLODE\_RK\_FWD (Radau1A LU)');
xlabel('Number of Steps');
ylabel('Relative Error');
title('Number of Steps vs. Relative Error (CBM4)');
saveas(gcf,strcat(resultsPath,'/CBM4_RK_ALL_Implicit_NstepsResults.eps'),'epsc');
saveas(gcf,strcat(resultsPath,'/CBM4_RK_ALL_Implicit_NstepsResults.fig'),'fig');
        