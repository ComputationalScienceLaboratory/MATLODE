delta = 0.000001;

Ode_Function = @flame_Function;
Ode_Jacobian = @flame_Jacobian;

Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'storeCheckpoint',true);

Tspan = [ 0; 2/delta ];
y0 = delta;

        % Compute Reference
        disp('Reference');
        options = odeset('AbsTol', 1d-13,'RelTol',1d-13,'Jacobian',Ode_Jacobian);
%         options = MATLODE_OPTIONS( Options, ...
%             'displaySteps', false, ...
%             'displayStats', false, ...
%             'storeCheckpoint', false, ...
%             'RelTol',     ones(2,1)*1d-13, ...
%             'AbsTol',     ones(2,1)*1d-13, ...
%             'Jacobian',   Ode_Jacobian, ...
%             'Max_no_steps', 400000, ...
%             'WarningConfig', 0);
       [ T_Ref, Y_Ref ] = ode15s(Ode_Function,Tspan,y0, options);
%         [ T_Ref, Y_Ref ] = ode113(Ode_Function,Tspan,y0, options);
%         [ T_Ref, Y_Ref ] = MATLODE_RK_FWD_Integrator(Ode_Function,Tspan,y0,options);
                
        Npoints = 11;
        TOLS = logspace(-11,-1,Npoints);
        
        ISTATUS_Array = cell(Npoints,1);
        RSTATUS_Array = cell(Npoints,1);

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
                              'RelTol',     ones(2,1)*RelTol, ...
                              'AbsTol',     ones(2,1)*AbsTol, ...
                              'Jacobian',   Ode_Jacobian, ...
                              'Max_no_steps', 400000, ...
                              'WarningConfig', 0);
            Options_MATLAB = odeset('AbsTol',AbsTol,'RelTol',RelTol,'Jacobian',Ode_Jacobian);
                          
%             tic;
%             output_ode45 = ode45(Ode_Function,Tspan,y0,Options_MATLAB);
%             elapsedTime_ode45(ipt) = toc;
%             T = transpose(output_ode45.x);
%             Y = transpose(output_ode45.y);
%             nsteps_ode45(ipt) = output_ode45.stats.nsteps;
%             error_ode45(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
            
            disp('ode15s');
            tic;
            output_ode15s = ode15s(Ode_Function,Tspan,y0,Options_MATLAB);
            elapsedTime_ode15s(ipt) = toc;
            T = transpose(output_ode15s.x);
            Y = transpose(output_ode15s.y);
            nsteps_ode15s(ipt) = output_ode15s.stats.nsteps;            
            error_ode15s(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
            
%             tic;
%             output_ode23 = ode23(Ode_Function,Tspan,y0,Options_MATLAB);
%             elapsedTime_ode23(ipt) = toc;
%             T = transpose(output_ode23.x);
%             Y = transpose(output_ode23.y);
%             nsteps_ode23(ipt) = output_ode23.stats.nsteps;            
%             error_ode23(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
%             
%             tic;
%             output_ode113 = ode113(Ode_Function,Tspan,y0,Options_MATLAB);
%             elapsedTime_ode113(ipt) = toc;
%             T = transpose(output_ode113.x);
%             Y = transpose(output_ode113.y);
%             nsteps_ode113(ipt) = output_ode113.stats.nsteps;            
%             error_ode113(ipt) = relativeError(Y(end,:),Y_Ref(end,:));            
            
            disp('ode23s');
            tic;
            output_ode23s = ode23s(Ode_Function,Tspan,y0,Options_MATLAB);
            elapsedTime_ode23s(ipt) = toc;
            T = transpose(output_ode23s.x);
            Y = transpose(output_ode23s.y);
            nsteps_ode23s(ipt) = output_ode23s.stats.nsteps;            
            error_ode23s(ipt) = relativeError(Y(end,:),Y_Ref(end,:));  
            
%             tic;
%             output_ode23t = ode23t(Ode_Function,Tspan,y0,Options_MATLAB);
%             elapsedTime_ode23t(ipt) = toc;
%             T = transpose(output_ode23t.x);
%             Y = transpose(output_ode23t.y);
%             nsteps_ode23t(ipt) = output_ode23t.stats.nsteps;            
%             error_ode23t(ipt) = relativeError(Y(end,:),Y_Ref(end,:));   
%             
%             tic;
%             output_ode23tb = ode23tb(Ode_Function,Tspan,y0,Options_MATLAB);
%             elapsedTime_ode23tb(ipt) = toc;
%             T = transpose(output_ode23tb.x);
%             Y = transpose(output_ode23tb.y);
%             nsteps_ode23tb(ipt) = output_ode23tb.stats.nsteps;            
%             error_ode23tb(ipt) = relativeError(Y(end,:),Y_Ref(end,:));             
            
%             disp('ERK_FWD');
%             Options = MATLODE_OPTIONS(Options,'Method',1);
%             [ T, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
%             elapsedTime_ERK_FWD_1(ipt) = RSTATUS.Etime;
%             nsteps_ERK_FWD_1(ipt) = ISTATUS.Nstp;
%             error_ERK_FWD_1(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
%             
%             disp('ERK_FWD');
%             Options = MATLODE_OPTIONS(Options,'Method',4);
%             [ T, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
%             elapsedTime_ERK_FWD_4(ipt) = RSTATUS.Etime;
%             nsteps_ERK_FWD_4(ipt) = ISTATUS.Nstp;
%             error_ERK_FWD_4(ipt) = relativeError(Y(end,:),Y_Ref(end,:));
%             
%             disp('ERK_FWD');
%             Options = MATLODE_OPTIONS(Options,'Method',5);
%             [ T, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
%             elapsedTime_ERK_FWD_5(ipt) = RSTATUS.Etime;
%             nsteps_ERK_FWD_5(ipt) = ISTATUS.Nstp;
%             error_ERK_FWD_5(ipt) = relativeError(Y(end,:),Y_Ref(end,:));            
            
            disp('ROS_FWD');
            Options = MATLODE_OPTIONS(Options,'Method',0);
            [ T, Y, Stats ] = MATLODE_ROS_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
            elapsedTime_ROS_FWD(ipt) = Stats.RSTATUS.Etime;
            nsteps_ROS_FWD(ipt) = Stats.ISTATUS.Nstp;
            error_ROS_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:));    
            
            disp('RK_FWD');
            [ T, Y, Stats ] = MATLODE_RK_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
            elapsedTime_RK_FWD(ipt) = Stats.RSTATUS.Etime;
            nsteps_RK_FWD(ipt) = Stats.ISTATUS.Nstp;
            error_RK_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:)); 
            
            disp('SDIRK_FWD');
            [ T, Y, Stats ] = MATLODE_SDIRK_FWD_Integrator( Ode_Function, Tspan, y0, Options );            
            elapsedTime_SDIRK_FWD(ipt) = Stats.RSTATUS.Etime;
            nsteps_SDIRK_FWD(ipt) = Stats.ISTATUS.Nstp;
            error_SDIRK_FWD(ipt) = relativeError(Y(end,:),Y_Ref(end,:));            
        end
        
        figure;
%         loglog(elapsedTime_ode45,error_ode45,'--o', ...
%                elapsedTime_ode15s,error_ode15s,'--+', ...
%                elapsedTime_ode23,error_ode23, '--*', ...
%                elapsedTime_ode113,error_ode113, '--^', ...
%                elapsedTime_ode23s,error_ode23s, '--x', ...
%                elapsedTime_ode23t,error_ode23t, '--s', ...
%                elapsedTime_ode23tb,error_ode23tb, '--d', ...
%                elapsedTime_ERK_FWD,error_ERK_FWD, '-o', ...
%                elapsedTime_ROS_FWD,error_ROS_FWD, '-+', ...
%                elapsedTime_RK_FWD,error_RK_FWD, '-*', ...
%                elapsedTime_SDIRK_FWD,error_SDIRK_FWD, '-^');
        loglog(elapsedTime_ode15s,error_ode15s,'--+', ...
               elapsedTime_ode23s,error_ode23s, '--x', ...
               elapsedTime_ROS_FWD,error_ROS_FWD, '-x', ...
               elapsedTime_RK_FWD,error_RK_FWD, '-s', ...
               elapsedTime_SDIRK_FWD,error_SDIRK_FWD, '-d');           
%                'ode15s', ...           
%                'ode23s', ...
%                'ode23t', ...
%                'ode23tb', ...                   
%                'MATLODE\_ROS\_FWD', ...
%                'MATLODE\_RK\_FWD', ...
%                'MATLODE\_SDIRK\_FWD');
%         legend('ode45', ...
%                'ode23', ...
%                'ode113', ...
%                'MATLODE\_ERK\_FWD (Erk23)', ...
%                'MATLODE\_ERK\_FWD (Dopri5)',...
%                'MATLODE\_ERK\_FWD (Dopri853)');
        legend('ode15s', ...
               'ode23s', ...
               'MATLODE\_ROS\_FWD (Lobatto3C)', ...
               'MATLODE\_RK\_FWD (Ros4)',...
               'MATLODE\_SDIRK\_FWD (Sdirk4A)');           
        xlabel('CPU Time (sec)');
        ylabel('Relative Error');
%         saveas(gcf,strcat(resultsPath,'/Lorenz3_Explicit_CpuResults.eps'),'epsc');
%         saveas(gcf,strcat(resultsPath,'/Lorenz3_Explicit_CpuResults.fig'),'fig');
         
        figure;
%         loglog(nsteps_ode45,error_ode45,'--o', ...
%                nsteps_ode15s,error_ode15s,'--+', ...
%                nsteps_ode23,error_ode23, '--*', ...
%                nsteps_ode113,error_ode113, '--^', ...
%                nsteps_ode23s,error_ode23s, '--x', ...
%                nsteps_ode23t,error_ode23t, '--s', ...
%                nsteps_ode23tb,error_ode23tb, '--d', ...
%                nsteps_ERK_FWD,error_ERK_FWD, '-o', ...
%                nsteps_ROS_FWD,error_ROS_FWD, '-+', ...
%                nsteps_RK_FWD,error_RK_FWD, '-*', ...
%                nsteps_SDIRK_FWD,error_SDIRK_FWD, '-^');
        loglog(nsteps_ode15s,error_ode15s,'--o', ...
            nsteps_ode23s,error_ode23s, '--x', ...
               nsteps_ROS_FWD,error_ROS_FWD, '-x', ...
               nsteps_RK_FWD,error_RK_FWD, '-s', ...
               nsteps_SDIRK_FWD,error_SDIRK_FWD, '-d');
%                'ode15s', ...           
%                'ode23s', ...
%                'ode23t', ...
%                'ode23tb', ...                   
%                'MATLODE\_ROS\_FWD', ...
%                'MATLODE\_RK\_FWD', ...
%                'MATLODE\_SDIRK\_FWD');
%         legend('ode45', ...
%                'ode23', ...
%                'ode15s', ...
%                'ode23s', ...
%                'ode113', ...
%                'MATLODE\_ERK\_FWD (Erk23)', ...
%                'MATLODE\_ERK\_FWD (Dopri5)',...
%                'MATLODE\_ERK\_FWD (Dopri853)');
        legend('ode15s', ...
               'ode23s', ...
               'MATLODE\_ROS\_FWD (Lobatto3C)', ...
               'MATLODE\_RK\_FWD (Ros4)',...
               'MATLODE\_SDIRK\_FWD (Sdirk4A)');  
        xlabel('Number of Steps');
        ylabel('Relative Error');
        title('Number of Steps vs. Relative Error (Van Der Pol)');
