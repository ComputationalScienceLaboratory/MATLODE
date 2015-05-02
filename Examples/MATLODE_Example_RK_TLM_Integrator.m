%% MATLODE_Example_RK_TLM_Integrator
%
% <html>
%    Up: <a href="../Examples/html/Examples.html">Examples</a>
% </html>
%%
% For the following examples Arenstorf Orbit is used as a toy problem 
% to illustrate |MATLODE_RK_TLM_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load the input parameters into the workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Ode_Jacobian        = @arenstorfOrbit_Jacobian;
Ode_YTLM            = eye(4);
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that the model is loaded in the workspace, one performs a
% forward implicit Runge-Kutta integration using the prebuilt default
% settings. Note, for implicite Runge-Kutta integrators, the Jacobian is
% required.
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Y_TLM',Ode_YTLM);
[ ~, Y, Sens ] = MATLODE_RK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Printing out our results, we can analyze our model state at our final
% time.
disp('solution at Time_Interval(2)');
disp(Y(end,:));
disp('sensitivity at Time_Interval(2)');
disp(Sens);

%% Advanced Features
% To perform *nondirect sensitivity analysis*, toggle the 'DirectTLM'
% option parameter to false. This enables the sensitivity matrix
% to be calculated using Newton iterations. Note, it is strongly
% recommended to first try direct sensitivity analysis before trying
% nondirect for efficiency purposes. 
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Y_TLM',Ode_YTLM,'DirectTLM',false);
[ ~, Y, Sens ] = MATLODE_RK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Printing out our results, we can analyze our model state at our final
% time.
disp('solution at Time_Interval(2)');
disp(Y(end,:));
disp('sensitivity at Time_Interval(2)');
disp(Sens);