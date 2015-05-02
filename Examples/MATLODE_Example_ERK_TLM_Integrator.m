%% MATLODE_Example_ERK_TLM_Integrator
%
% <html>
%    Up: <a href="../Examples/html/Examples.html">Examples</a>
% </html>
%%
% For the following examples we will use Arenstorf Orbit as a toy problem 
% to illustrate |MATLODE_ERK_TLM_Integrator| functionalities and features.
% To initially setup Brusselator, execute the MATLAB commands below to
% load our input parameters into our workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Ode_Jacobian        = @arenstorfOrbit_Jacobian;
Ode_YTLM            = eye(4);
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that we have our model loaded in our workspace, we can perform a
% tangent linear explicit Runge-Kutta integration using MATLODE's prebuilt default
% settings. We note that a Jacobian and Y_TLM required and passed by MATLODE's
% option struct.
Options  = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Y_TLM',Ode_YTLM);
[ ~, Y, Sens ] = MATLODE_ERK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Printing out our results, we can analyze our model state at our final
% time.
disp('solution at Time_Interval(2)');
disp(Y(end,:));
disp('sensitivity at Time_Interval(2)');
disp(Sens);

%% Advanced Features
% To *save the model state at each time step*, one needs to
% initialize a MATLODE(R) option struct to store the fine tuning settings. The
% |(key,value)| pair associated for saving the model state at each time
% step is denoted as |('storeCheckpoint',true)| or
% |('storeCheckpoint',false)| depending on whether or not one wants to
% explicitly fine tune the integrator. In this case, the intermediary time
% step values are stored executing the command below.
Options = MATLODE_OPTIONS('storeCheckpoint',true,'Jacobian',Ode_Jacobian,'Y_TLM',Ode_YTLM);

%%
% To run |MATLODE_ERK_FWD_Integrator| using the fine tuning, one needs to
% insert the option struct into the integrator's fourth parameter position.
[ ~, Y, Sens ] = MATLODE_ERK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Printing out our results, we can analyze our model state at our final
% time.
format long;
disp('solution at Time_Interval(2)');
disp(Y(end,:));
disp('sensitivity at Time_Interval(2)');
disp(Sens);

%%
% After plotting the results, one can now visualize the model.
figure(1);
plot(Y(:,1),Y(:,2));
title('Arenstorf Orbit');
xlabel('Y(:,1)');
ylabel('Y(:,2)');

%%
% To obtain a smoother graphical repusentation, one can further tighten the 
% error tolerances. To tighten the relative and abolute error tolerances,
% one fine tunes the option struct. Since the option struct is already in 
% the workspace, one adds the relative and absolute |(key,value)| pair to 
% the option struct. Then plot the results.
Options = MATLODE_OPTIONS(Options,'AbsTol',1e-12,'RelTol',1e-12);
[ T, Y, Sens ] = MATLODE_ERK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Printing out our results, we can analyze our model state at our final
% time.
format long;
disp('solution at Time_Interval(2)');
disp(Y(end,:));
disp('sensitivity at Time_Interval(2)');
disp(Sens);

%%
% After plotting the results, one can now visualize the model.
figure(2);
plot(Y(:,1),Y(:,2));
title('Arenstorf Orbit');
xlabel('Y(:,1)');
ylabel('Y(:,2)');
