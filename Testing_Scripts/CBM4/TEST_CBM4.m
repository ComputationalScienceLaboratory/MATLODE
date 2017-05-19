clear global;

% analysis: -1: solution
%           -2: error
%           -3: cpu time
analysis = -3;

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

switch ( analysis )
    case -1 % solution
        switch ( integrator )
            case -1
                disp( 'Solving problem with MATLODE_ERK_FWD_Integrator: ' );
                [ T, Y, Stats] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );
                
                implementation = 'FWD';
                family = 'ERK';
                           
            case -2
                disp( 'Solving problem with MATLODE_ROS_FWD_Integrator: ' );
                [ T, Y, Stats ] = MATLODE_ROS_FWD_Integrator( Ode_Function, Tspan, y0, Options  );  
                
                implementation = 'FWD';
                family = 'ROS';                 

            case -3
                disp( 'Solving problem with MATLODE_RK_FWD_Integrator: ' );
                [ T, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );  
                
                implementation = 'FWD';
                family = 'RK';                  

            case -4
                disp( 'Solving problem with MATLODE_SDIRK_FWD_Integrator: ');
                [ T, Y, Stats ] = MATLODE_SDIRK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );
                
                implementation = 'FWD';
                family = 'SDIRK';                  

            case -5
                disp( 'Solving problem with MATLODE_ERK_TLM_Integrator: ' );
                [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_ERK_TLM_Integrator( Ode_Function, Tspan, y0, Options );
  
                implementation = 'TLM';
                family = 'ERK';
                
            case -6 
                disp( 'Solving problem with ROS_TLM_Integrator: ');
                [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_ROS_TLM_Integrator( Ode_Function, Tspan, y0, Options );               
                
                implementation = 'TLM';
                family = 'ROS';
                
            case -7
                disp( 'Solving problem with RK_TLM_Integrator: ');
                [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_RK_TLM_Integrator( Ode_Function, Tspan, y0, Options );
                
                implementation = 'TLM';
                family = 'RK';
 
            case -8
                disp( 'Solving problem with SDIRK_TLM_Integrator: ');
                [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_SDIRK_TLM_Integrator( Ode_Function, Tspan, y0, Options );
                
                implementation = 'TLM';
                family = 'SDIRK';

            case -9
                disp( 'Solving problem with ERK_ADJ_Integrator: ');
                [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_ERK_ADJ_Integrator( Ode_Function, Tspan, y0, Options );
                
                implementation = 'ADJ';
                family = 'ERK';                

            case -10
                disp( 'Solving problem with ROS_ADJ_Integrator: ');
                [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_ROS_ADJ_Integrator( Ode_Function, Tspan, y0, Options );
                
                implementation = 'ADJ';
                family = 'ROS';

            case -11
                disp( 'Solving problem with RK_ADJ_Integrator: ');
                [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_RK_ADJ_Integrator( Ode_Function, Tspan, y0, Options);
                
                implementation = 'ADJ';
                family = 'RK';                

            case -12
                disp( 'Solving problem with SDIRK_ADJ_Integrator: ');
                [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                    MATLODE_SDIRK_ADJ_Integrator( Ode_Function, Tspan, y0, Options );
                
                implementation = 'ADJ';
                family = 'SDIRK';   

        end
        
        for i = 1:NMONITOR
            figure; plot( (T)/3600, Y(:,MONITOR(i))/CFACTOR );
            title( SPC_NAMES( MONITOR(i),:) ,'FontSize',12); 
            set(gca,'XLim',[TSTART,TEND]/3600,'FontSize',12);
            xlabel('Time [ h ]','FontSize',12);
            ylabel('Concentration','FontSize',12);
        end  
        
    case -2 % analysis
        Npoints = 11;
        err = zeros(Npoints,1);
        steps = zeros(Npoints,1);
        TOLS = logspace(-11,-1,Npoints);
        
        ISTATUS_Array = cell(Npoints,1);
        RSTATUS_Array = cell(Npoints,1);

        for ipt=1:Npoints
            RelTol=TOLS(ipt);
            AbsTol = RelTol;
            Options = MATLODE_OPTIONS( Options, ...
                              'displaySteps', 0, ...
                              'storeCheckpoint', 0, ...
                              'RelTol',     ones(yDimension,1)*RelTol, ...
                              'AbsTol',     ones(yDimension,1)*AbsTol, ...
                              'Jacobian',   Ode_Jacobian, ...
                              'AbsTol_TLM', ones(1,yDimension)*AbsTol, ...
                              'RelTol_TLM', ones(1,yDimension)*RelTol, ...
                              'Y_TLM',      eye(yDimension,yDimension), ...
                              'NTLM',       yDimension, ...
                              'AbsTol_ADJ', ones(1,yDimension)*AbsTol, ...
                              'AbsTol_ADJ', ones(1,yDimension)*RelTol, ...
                              'WarningConfig', 0, ...
                              'NADJ',       yDimension );
            switch ( integrator )
                case -1
                    disp( 'Solving problem with MATLODE_ERK_FWD_Integrator: ' );
                    [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );

                    implementation = 'FWD';
                    family = 'ERK';

                case -2
                    disp( 'Solving problem with MATLODE_ROS_FWD_Integrator: ' );
                    [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ROS_FWD_Integrator( Ode_Function, Tspan, y0, Options  );  

                    implementation = 'FWD';
                    family = 'ROS';                 

                case -3
                    disp( 'Solving problem with MATLODE_RK_FWD_Integrator: ' );
                    [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );  

                    implementation = 'FWD';
                    family = 'RK';                  

                case -4
                    disp( 'Solving problem with MATLODE_SDIRK_FWD_Integrator: ');
                    [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_SDIRK_FWD_Integrator( Ode_Function, Tspan, y0, Options  );

                    implementation = 'FWD';
                    family = 'SDIRK';                  

                case -5
                    disp( 'Solving problem with MATLODE_ERK_TLM_Integrator: ' );
                    [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ERK_TLM_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'TLM';
                    family = 'ERK';

                case -6 
                    disp( 'Solving problem with ROS_TLM_Integrator: ');
                    [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ROS_TLM_Integrator( Ode_Function, Tspan, y0, Options );               

                    implementation = 'TLM';
                    family = 'ROS';

                case -7
                    disp( 'Solving problem with RK_TLM_Integrator: ');
                    [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_RK_TLM_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'TLM';
                    family = 'RK';

                case -8
                    disp( 'Solving problem with SDIRK_TLM_Integrator: ');
                    [ T, Y, Y_TLM, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_SDIRK_TLM_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'TLM';
                    family = 'SDIRK';

                case -9
                    disp( 'Solving problem with ERK_ADJ_Integrator: ');
                    [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ERK_ADJ_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'ADJ';
                    family = 'ERK';                

                case -10
                    disp( 'Solving problem with ROS_ADJ_Integrator: ');
                    [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_ROS_ADJ_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'ADJ';
                    family = 'ROS';

                case -11
                    disp( 'Solving problem with RK_ADJ_Integrator: ');
                    [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_RK_ADJ_Integrator( Ode_Function, Tspan, y0, Options);

                    implementation = 'ADJ';
                    family = 'RK';                

                case -12
                    disp( 'Solving problem with SDIRK_ADJ_Integrator: ');
                    [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                        MATLODE_SDIRK_ADJ_Integrator( Ode_Function, Tspan, y0, Options );

                    implementation = 'ADJ';
                    family = 'SDIRK';  

            end
            T_Array{ipt} = T;
            Y_Array{ipt} = Y;
            ISTATUS_Array{ipt} = ISTATUS;
            RSTATUS_Array{ipt} = RSTATUS;
            Ierr_Array{ipt} = Ierr;
            
            disp('Solving problem with ode15s: ')
            o = odeset('AbsTol',1d-12,'RelTol',1d-12);
            if ( ipt == 1 )
                [t,y] = ode45(Ode_Function,Tspan,y0, o);
            end
            switch ( error_mode )
                case -1
                    disp('relative root mean squared (rRMS)');
                    err_matlab(ipt,1) = rRMS(Y(length(T),:),y(length(t),:));
                case -2
                    disp('relative error');
                    err_matlab(ipt,1) = relativeError(Y(length(T),:),y(length(t),:));
                case -3
                    disp('absolute error');
                    err_matlab(ipt,1) = absoluteError(Y(length(T),:),y(length(t),:));
                otherwise
                    disp('Error: Choose error_mode');
            end
            
            steps(ipt) = ISTATUS.Nstp;
        end
        
        saveVariable( [resultsPath, '/CBM4_Results_', implementation, ...
            '/CBM4_Results_', implementation, '_', family, '/', Coefficient.Name, '/'], ...
            ['Error_' family '_' implementation '_CBM4_T'], T_Array, ...
            ['Error_' family '_' implementation '_CBM4_Y'], Y_Array, ...
            ['Error_' family '_' implementation '_CBM4_ISTATUS'], ISTATUS_Array, ...
            ['Error_' family '_' implementation '_CBM4_RSTATUS'], RSTATUS_Array, ...
            ['Error_' family '_' implementation '_CBM4_Ierr'], Ierr_Array, enableSave );
        
        figure;
        loglog(steps,err_matlab);
        [convergentRate, ~ ] = polyfit(log(steps(1:3)),log(err_matlab(1:3)),1);
        titleName = {['Error: ' family ' ' implementation ' Integrator (CBM4)'], ...
                ['Coefficient: ' num2str(Coefficient.Name)], ...
                ['Rate of Convergent: ', num2str(convergentRate(1))] };
        title(titleName);
        ylabel('Relative Error');
        xlabel('Number of steps');
        if ( enableSave == true )
            saveas(gcf, [resultsPath, '/CBM4_Results_', implementation, ...
                '/CBM4_Results_', implementation, '_', family, '/', Coefficient.Name, ...
                '/Error_', family, '_', implementation,'_CBM4_', Coefficient.Name, '.fig']);       
        end
        
case -3
        % Compute Reference
        disp('Reference ode15s');
        options = odeset('AbsTol', 1d-8,'RelTol',1d-8,'Jacobian',Ode_Jacobian);
       [ ~, Y_Ref ] = ode15s(Ode_Function,Tspan,y0, options);
       
        disp('Reference ode15s');
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
        legend('ode15s', ...          
               'ode23s', ...
               'ode23t', ...
               'ode23tb', ...
               rosTitle, ...
               rkTitle, ...
               sdirkTitle);             
        xlabel('CPU Time (sec)');
        ylabel('Relative Error');
        title('CPU vs. Relative Error (CBM4)');
        saveas(gcf,strcat(resultsPath,'/CBM4_Implicit_CpuResults.eps'),'epsc');
        saveas(gcf,strcat(resultsPath,'/CBM4_Implicit_CpuResults.fig'),'fig');
         
        figure;
        loglog(nsteps_ode15s,error_ode15s,'--+', ...
               nsteps_ode23s,error_ode23s, '--x', ...
               nsteps_ode23t,error_ode23t, '--s', ...
               nsteps_ode23tb,error_ode23tb, '--d', ...
               nsteps_ROS_FWD,error_ROS_FWD, '-+', ...
               nsteps_RK_FWD,error_RK_FWD, '-*', ...
               nsteps_SDIRK_FWD,error_SDIRK_FWD, '-^');        
        legend('ode15s', ...          
               'ode23s', ...
               'ode23t', ...
               'ode23tb', ...
               rosTitle', ...
               rkTitle, ...
               sdirkTitle);     
        xlabel('Number of Steps');
        ylabel('Relative Error');
        title('Number of Steps vs. Relative Error (CBM4)');
        saveas(gcf,strcat(resultsPath,'/CBM4_Implicit_NstepsResults.eps'),'epsc');           
        saveas(gcf,strcat(resultsPath,'/CBM4_Implicit_NstepsResults.fig'),'fig');
        
    otherwise
        disp( 'Error: Choose analysis' );
end

% profreport('TEST_CBM4')