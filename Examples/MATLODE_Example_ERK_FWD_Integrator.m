%% MATLODE_Example_ERK_FWD_Integrator
%
% <html>
%    Up: <a href="../Examples/html/Examples.html">Examples</a>
% </html>
%%
% For the following examples Arenstorf Orbit is used as a toy problem 
% to illustrate |MATLODE_ERK_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load the input parameters into the workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that the model is loaded in the workspace, one performs a
% forward explicit Runge-Kutta integration using the prebuilt default
% settings.
[ ~, Y ] = MATLODE_ERK_FWD_Integrator(Ode_Function,Time_Interval,Y0);

%%
% Execute the following commands to analyze the final model state.
disp('solution at Time_Interval(2)');
disp(Y(end,:));

%% Advanced Features
% To *save the model state at each time step*, one needs to
% initialize a MATLODE(R) option struct to store the fine tuning settings. The
% |(key,value)| pair associated for saving the model state at each time
% step is denoted as |('storeCheckpoint',true)| or
% |('storeCheckpoint',false)| depending on whether or not one wants to
% explicitly fine tune the integrator. In this case, the intermediary time
% step values are stored executing the command below.
Options = MATLODE_OPTIONS('storeCheckpoint',true);

%%
% To run |MATLODE_ERK_FWD_Integrator| using the fine tuning, one needs to
% insert the option struct into the integrator's fourth parameter position.
[ ~, Y ] = MATLODE_ERK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

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
[ T, Y ] = MATLODE_ERK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);
figure(2);
plot(Y(:,1),Y(:,2));
title('Arenstorf Orbit');
xlabel('Y(:,1)');
ylabel('Y(:,2)');


