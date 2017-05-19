% Testing Suite: Van Der Pol

%clear global;

% integrator: -1: ERK_FWD_DRIVER_Integrator
%             -2: ROS_FWD_DRIVER_Integrator
%             -3: RK_FWD_DRIVER_Integrator
%             -4: SDIRK_FWD_DRIVER_Integrator
%
%             -5: ERK_TLM_DRIVER_Integrator
%             -6: ROS_TLM_DRIVER_Integrator
%             -7: RK_TLM_DRIVER_Integrator
%             -8: SDIRK_TLM_DRIVER_Integrator
%
%              -9: ERK_ADJ_DRIVER_Integrator
%             -10: ROS_ADJ_DRIVER_Integrator
%             -11: RK_ADJ_DRIVER_Integrator
%             -12: SDIRK_ADJ_DRIVER_Integrator              
integrator = -1;

% Adjoint: -1: xxxxADJ1 w/ Quadrature
%          -2: xxxxADJ1 no Quadrature
%          -3: xxxxADJ2 w/ Quadrature
%          -4: xxxxADJ2 no Quadrature
adj_mode = -1;

disp( 'Ver Der Pol Problem: ' );
Tspan = [0 500];
y0 = [2;-0.66];

% SDIRK Coefficient     ERK Coefficient     RK Coefficient      Ros Coefficient
%   1: Sdirk2A              1: Erk23            1: Rada2A           1: Ros2
%   2: Sdirk2B              2: Erk3_Heun        2: Lobatto3C        2: Ros3
%   3: Sdirk3A              3: Erk43            3: Gauss            3: Ros4
%   4: Sdirk4A              4: Dopri5           4: Radau1A          4: Rodas3
%   5: Sdirk4B              5: Verme            5: Lobatto3A        5: Rodas4
%                           6: Dopri853

RelTol = [ 1e-5 1e-5 ];
AbsTol = [ 1e-5 1e-5 ];

ICNTRL = fatOde_ICNTRL_Set( 'Autonomous',   0, ... % 1
                            'ITOL',         0, ... % 2
                            'Method',       0, ... % 3
                            'Max_no_steps', 0, ... % 4
                            'NewtonMaxit',  0, ... % 5
                            'StartNewton',  0, ... % 6
                            'DirectTLM',    0, ... % 7
                            'SaveLU',       0, ... % 8                        
                            'TLMNewtonEst', 0, ... % 9
                            'SdirkError',   0, ... % 10
                            'Gustafsson',   0, ... % 11 % 1: fwd_rk
                            'TLMtruncErr',  0, ... % 12
                            'FDAprox',      0, ... % 13                       
                            'AdjointSolve', 0, ... % 14
                            'AdjointType',  0, ... % 15
                            'DirectADJ',    0, ... % 16
                            'ChunkSize',    50 );  % 17

RCNTRL = fatOde_RCNTRL_Set( 'Hmin',      0, ...
                            'Hmax',      0, ...  
                            'Hstart',    0, ...
                            'FacMin',    0, ...
                            'FacMax',    0, ...
                            'FacRej',    0, ...
                            'FacSafe',   0, ...
                            'ThetaMin',  0, ...
                            'NewtonTol', 0, ...
                            'Qmin',      0, ...
                            'Qmax',      0 );

Options = fatOde_OPTIONS_Set( 'AbsTol',       AbsTol, ...
                              'RelTol',       RelTol, ...
                              'Jacobian',     @vanDerPol_Jacobian, ...
                              'AbsTol_TLM',   AbsTol.*10, ...
                              'RelTol_TLM',   RelTol.*10, ...
                              'Y_TLM',        eye(2,2), ...
                              'NTLM',         2, ...
                              'AbsTol_ADJ',   AbsTol.*10, ...
                              'RelTol_ADJ',   RelTol.*10, ...
                              'NADJ',         2, ...
                              'Desired_Mode', -1*adj_mode, ...
                              'displayStats', true, ...
                              'displaySteps', false );

Npoints = 7;
err = zeros(Npoints,1);
steps = zeros(Npoints,1);
TOLS = logspace(-7,-1,Npoints);    
for ipt=1:Npoints
    RelTol=TOLS(ipt);
    AbsTol = RelTol;
    Options = fatOde_OPTIONS_Set( 'RelTol',     ones(1,2)*RelTol, ...
                      'AbsTol',     ones(1,2)*AbsTol, ...
                      'Jacobian',   @vanDerPol_Jacobian, ...
                      'AbsTol_TLM', ones(1,2)*AbsTol, ...
                      'RelTol_TLM', ones(1,2)*RelTol, ...
                      'Y_TLM',      eye(2,2), ...
                      'NTLM',       2, ...
                      'AbsTol_ADJ', ones(1,2)*AbsTol, ...
                      'AbsTol_ADJ', ones(1,2)*RelTol, ...
                      'NADJ',       2 );
    switch ( integrator )
        case -1
            disp( 'Solving problem with ERK_FWD_DRIVER_Integrator: ' );
%profile clear;
%profile on;     
            cput = cputime;
            [ T_1, Y_1, ISTATUS_1, RSTATUS_1, Ierr_1, Coefficient_1 ] = ...
                ERK_FWD_DRIVER_Integrator_old( @vanDerPol_Function, Tspan, y0, Options, RCNTRL, ICNTRL );
            cputime_matlODE_old(ipt) = cputime - cput;  
%            cput = cputime;              
%            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
%                ERK_FWD_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, Options, RCNTRL, ICNTRL );
%            cputime_matlODE(ipt) = cputime - cput;            
       
%profreport('TEST_vanDerPol')
            implementation = 'FWD';
            family = 'ERK';

%             Options.AbsTol = Options.AbsTol(1);
%             Options.RelTol = Options.RelTol(1);
%             cput = cputime;
%             [t,y] = ode45(@vanDerPol_Function,Tspan,y0, Options);  
%             cputime_matlab(ipt) = cputime - cput;

        case -2
            disp( 'Solving problem with ROS_FWD_DRIVER_Integrator: ' );
            cput = cputime;
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ROS_FWD_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );   
            cputime_matlODE(ipt) = cputime - cput;

            implementation = 'FWD';
            family = 'ROS';

            Options.AbsTol = Options.AbsTol(1);
            Options.RelTol = Options.RelTol(1);
            cput = cputime;
            [t,y] = ode23(@vanDerPol_Function,Tspan,y0, Options);                  
            cputime_matlab(ipt) = cputime - cput;

        case -3
            disp( 'Solving problem with RK_FWD_DRIVER_Integrator: ' );
            cput = cputime;
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = RK_FWD_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );
            cputime_matlODE(ipt) = cputime - cput;

            implementation = 'FWD';
            family = 'RK';

            Options.AbsTol = Options.AbsTol(1);
            Options.RelTol = Options.RelTol(1);
            cput = cputime;
            [t,y] = ode15s(@vanDerPol_Function,Tspan,y0, Options);                  
            cputime_matlab(ipt) = cputime - cput;

        case -4
            disp( 'Solving problem with SDIRK_FWD_DRIVER_Integrator: ');
            cput = cputime;
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = SDIRK_FWD_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );
            cputime_matlODE(ipt) = cputime - cput;

            implementation = 'FWD';
            family = 'SDIRK';

            Options.AbsTol = Options.AbsTol(1);
            Options.RelTol = Options.RelTol(1);
            cput = cputime;
            [t,y] = ode15s(@vanDerPol_Function,Tspan,y0, Options);                  
            cputime_matlab(ipt) = cputime - cput;

        case -5
            disp( 'Solving problem with ERK_TLM_DRIVER_Integrator: ' );
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ERK_TLM_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'TLM';
            family = 'ERK';
        case -6 
            disp( 'Solving problem with ROS_TLM_DRIVER_Integrator: ');
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                ROS_TLM_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, Options, RCNTRL, ICNTRL );               

            implementation = 'TLM';
            family = 'ROS';

        case -7
            disp( 'Solving problem with RK_TLM_DRIVER_Integrator: ');
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = ...
                RK_TLM_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, Options, RCNTRL, ICNTRL );

            implementation = 'TLM';
            family = 'RK';

        case -8
            disp( 'Solving problem with SDIRK_TLM_DRIVER_Integrator: ');
            [ T, Y, ISTATUS, RSTATUS, Ierr, Coefficient ] = SDIRK_TLM_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'TLM';
            family = 'SDIRK';

        case -9
            disp( 'Solving problem with ERK_ADJ_DRIVER_Integrator: ');
            [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ERK_ADJ_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'ADJ';
            family = 'ERK';                

        case -10
            disp( 'Solving problem with ROS_ADJ_DRIVER_Integrator: ');
            [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = ROS_ADJ_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'ADJ';
            family = 'ROS';

        case -11
            disp( 'Solving problem with RK_ADJ_DRIVER_Integrator: ');
            [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = RK_ADJ_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'ADJ';
            family = 'RK';                

        case -12
            disp( 'Solving problem with SDIRK_ADJ_DRIVER_Integrator: ');
            [ T, Y, Lambda, ISTATUS, RSTATUS, Ierr, Coefficient ] = SDIRK_ADJ_DRIVER_Integrator( @vanDerPol_Function, Tspan, y0, ...
                Options, RCNTRL, ICNTRL );

            implementation = 'ADJ';
            family = 'SDIRK';   

    end

%     if ( Ierr(:) >= 0.0 )
%         titleName = {['Solution: ' family ' ' implementation ' Integrator (Van Der Pol)'], ...
%             ['Coefficient: ' num2str(Coefficient.Name) ] };        
% 
%         plot( T, Y );
%         xlabel('Time');
%         ylabel('Y1 & Y2');
%         title(titleName);
%     end   
%     steps(ipt) = length(T);    
end
semilogy(cputime_matlODE,TOLS,cputime_matlODE_old,TOLS);
title('Vectorization');
ylabel('ROLS');
xlabel('CPU Time');
legend('optimized','old version');