%% MATLODE_Example_ROS_FWD_Integrator
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
% to illustrate |MATLODE_ROS_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load the input parameters into the workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Ode_Jacobian        = @arenstorfOrbit_Jacobian;
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that the model is loaded in the workspace, one performs a
% forward Rosenbrock integration using the prebuilt default
% settings. Note, for Rosenbrock integrators, the Jacobian is
% required.
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian);
[ ~, Y ] = MATLODE_ROS_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Execute the following commands to analyze the final model state.
disp('solution at Time_Interval(2)');
disp(Y(end,:));

%% Advanced Features
% To analyze intermediary steps, toggle the option struct parameter
% 'displaySteps' to true. Printing intermediary steps often gives an
% immediate visual queue on how hard the intergrator is working internally.
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'displaySteps',true);
[ ~, Y ] = MATLODE_ROS_FWD_Integrator(Ode_Function,[Time_Interval(1) Time_Interval(1)+0.0029],Y0,Options);

%%
% Noticing above that no initial steps are rejected, increasing Hstart will
% force the error controller to be more aggressive in the beginning of the
% Rosenbrock intergration scheme. 
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'displaySteps',true,'Hstart',0.0005);
[ ~, Y ] = MATLODE_ROS_FWD_Integrator(Ode_Function,[Time_Interval(1) Time_Interval(1)+0.0029],Y0,Options);

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