%% MATLODE_Example_RK_FWD_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
% <html>
%    Up: <a href="../Examples/html/Examples.html">Examples</a>
% </html>
%%
% For the following examples Arenstorf Orbit is used as a toy problem 
% to illustrate |MATLODE_RK_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load the input parameters into the workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Ode_Jacobian        = @arenstorfOrbit_Jacobian;
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that the model is loaded in the workspace, one performs a
% forward implicit Runge-Kutta integration using the prebuilt default
% settings. Note, for implicite Runge-Kutta integrators, the Jacobian is
% required.
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian);
[ ~, Y ] = MATLODE_RK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

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
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'storeCheckpoint',true);
[ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% After plotting the results, one can now visualize the model.
figure(1);
plot(Y(:,1),Y(:,2));
title('Arenstorf Orbit');
xlabel('Y(:,1)');
ylabel('Y(:,2)');
PrintISTATUS(Stats.ISTATUS);

%%
% Depending on the model, it may be advantageous to use Kjell Gustafsson's
% error control approach described in [1]. One notes that in this toy
% problem, it is better to use Gusafsson, but for illustration purposes
% Gustafsson is toggled off to demonstrate the effect. 
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'storeCheckpoint',true,'Gustafsson',false);
[ ~, Y, Stats ] = MATLODE_RK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% After plotting the results and printing the integrator statistics, one 
% can visualize the model and compare statistics to the previous example.
figure(2);
plot(Y(:,1),Y(:,2));
title('Arenstorf Orbit');
xlabel('Y(:,1)');
ylabel('Y(:,2)');
PrintISTATUS(Stats.ISTATUS);


%% Reference
% [1] K. Gustafsson,  Control of error and convergence in ODE solvers,  1992 :Dept. of Automat. Contr., Lund Inst. Technol.

%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%

%%
% <html>
%   <div>
%       <img style="float: right" src="../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html> 