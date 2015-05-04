%% MATLODE_Example_SDIRK_FWD_Integrator
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
% to illustrate |MATLODE_SDIRK_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load the input parameters into the workspace. 
Ode_Function        = @arenstorfOrbit_Function;
Ode_Jacobian        = @arenstorfOrbit_Jacobian;
Ode_JacobianVector  = @arenstorfOrbit_JacobianVector;
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

%% Basic Functionality
% Now that the model is loaded in the workspace, one performs a
% forward Rosenbrock integration using the prebuilt default
% settings. 
Options = MATLODE_OPTIONS('Jacobian',Ode_Jacobian);
[ ~, Y ] = MATLODE_SDIRK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Execute the following commands to analyze the final model state.
disp('solution at Time_Interval(2)');
disp(Y(end,:));

%% Advanced Functionality
% If an analytical Jacobian is not available, MATLODE(R) can approximate
% the Jacobian by toggling the 'MatrixFree' option parameter. 
Options = MATLODE_OPTIONS('MatrixFree',true);
[ ~, Y ] = MATLODE_SDIRK_FWD_Integrator(Ode_Function,[Time_Interval(1) Time_Interval(1)+0.1],Y0,Options);

%%
% Execute the following commands to analyze the final model state.
disp('solution at Time_Interval(2)');
disp(Y(end,:));

%%
% If a Jacobian vector product approximation is available, can pass the
% Jacobian vector product function handler to 'Jacobian' and toggle the
% 'MatrixFree' in the option struct. 
Options = MATLODE_OPTIONS('Jacobian',Ode_JacobianVector,'MatrixFree',true);
[ ~, Y ] = MATLODE_SDIRK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);

%%
% Execute the following commands to analyze the final model state.
disp('solution at Time_Interval(2)');
disp(Y(end,:));

%%
% Copyright 2015 Computational Science Laboratory

%%
% <html>
%   <div>
%       <img style="float: right" src="../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>