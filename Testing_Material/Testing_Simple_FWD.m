%% Testing Material: Testing_Simple_FWD
%
% <html>
% Up: <a href="Testing_Material.html">Testing Material</a>
% </html>
% 
%%
% The forward test cases make ose of the following ODE test problem. 
Ode_Function = @vanDerPol_Function;
Ode_Jacobian = @vanDerPol_Jacobian;

Tspan = [0;  20  ];
y0    = [2; -0.66];

%%
% We first generate a reference solution using MATLAB's built in ODE
% integrator ode45.
Option_Ref = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[T_Ref, Y_Ref] = ode45(Ode_Function,Tspan,y0,Option_Ref);

%% Methodology
% In order not overwhelm the reader with very many test cases, we only
% demonstrate selecting the most common parameters in MATLODE(R). For all 
% ODE solvers in MATLODE(R) forward package, we perform an absolute and 
% relative tolerance sweep using the built in default ODE option struct 
% denoted by 'A' in the section title. Additionaly, we perform a second 
% absolute and relative tolerance sweep, denoted by 'B', changing the 
% method used in the integrator. In conclusion, this reduces the amount of 
% test cases for the reader substantially while illustrating the toolbox's 
% expective behavior. 

%% Test Case 1A: MATLODE_ERK_FWD_Integrator (Default)
% *Description*: Running |MATLODE_ERK_FWD_Integrator| with different absolute
% and relative tolerance values given default option struct.
Option = MATLODE_OPTIONS();
[steps, errorSolution] = TEST_General_Forward(@MATLODE_ERK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 1B: MATLODE_ERK_FWD_Integrator
% *Description*: Running |MATLODE_ERK_FWD_Integrator| with different absolute
% and relative tolerance values given alternative coefficients.
Option = MATLODE_OPTIONS('Method',3);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_ERK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 2A: MATLODE_RK_FWD_Integrator (Default)
% *Description*: Running |MATLODE_RK_FWD_Integrator| with different absolute
% and relative tolerance values given default option struct.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_RK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 2B: MATLODE_RK_FWD_Integrator 
% *Description*: Running |MATLODE_RK_FWD_Integrator| with different absolute
% and relative tolerance values given alternative coefficients.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Method',1);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_RK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 3A: MATLODE_ROS_FWD_Integrator (Default)
% *Description*: Running |MATLODE_ROS_FWD_Integrator| with different absolute
% and relative tolerance values given default option struct.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Max_no_steps',400000);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_ROS_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 3B: MATLODE_ROS_FWD_Integrator
% *Description*: Running |MATLODE_ROS_FWD_Integrator| with different absolute
% and relative tolerance values given alternative coefficients.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Max_no_steps',400000,'Method',2);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_ROS_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 4A: MATLODE_SDIRK_FWD_Integrator (Default)
% *Description*: Running |MATLODE_SDIRK_FWD_Integrator| with different absolute
% and relative tolerance values given default option struct.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_SDIRK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');

%% Test Case 4B: MATLODE_SDIRK_FWD_Integrator
% *Description*: Running |MATLODE_SDIRK_FWD_Integrator| with different absolute
% and relative tolerance values given alternative coefficients.
Option = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Method',5);
[steps, errorSolution] = TEST_General_Forward(@MATLODE_SDIRK_FWD_Integrator,Ode_Function,Tspan,y0,Option,-11,-4,T_Ref,Y_Ref);
loglog(steps,errorSolution);
title('Van Der Pol (Solution)'); ylabel('RMS Relative Error'); xlabel('Number of Steps');