%% MATLODE_SDIRK_TLM_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%                                 MATLODE_SDIRK_TLM_Integrator
%                [ T, Y, Sens ] = MATLODE_SDIRK_TLM_Integrator(Ode_Function, Time_Interval, Y0, Options)
%   [ T, Y, Sens, Quad, Stats ] = MATLODE_SDIRK_TLM_Integrator(Ode_Function, Time_Interval, Y0, Options)
%
%% Input Parameters
% |Ode_Function|: model function
%
% |Time_Interval|: time span
%
% |Y0|: initial model state vector
%
% |Options|: MATLODE option struct
%
%% Output Parameters
% |T|: saved time snapshots
%
% |Y|: saved model state vectors
%
% |Sens|: Sensitivity matrix
%
% |Quad|: Quadrature term
%
% |Stats|: integrator statistics
%
%% Description
% Driver file to solve the system y' = F(t,y) and tangent linear sensitivity 
% using an Singly Diagonally Implicit Runge-Kutta (SDIRK) method.
%
% |MATLODE_SDIRK_TLM_Integrator| displays the available methods
% associated with the exponential forward integrator.
%
% |[T, Y, Sens] = MATLODE_SDIRK_TLM_Integrator(Ode_Function, Time_Interval, Y0,
% Options)| computes the ODE solution with respect to the user supplied
% options configuration and tangent linear sensitivity.
%
% |[T, Y, Sens, Quad, Stats]| = MATLODE_SDIRK_TLM_Integrator(Ode_Function,
% Time_Interval, Y0, Options)| computes the ODE solution with respect to 
% the user supplied options configuration, tangent linear sensitivity,
% quadrature and statisitics.
%
%% Example
% For the following examples we will use Arenstorf Orbit as a toy problem 
% to illustrate |MATLODE_SDIRK_TLM_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load our input parameters into our workspace. 
%
%   Ode_Function        = @arenstorfOrbit_Function;
%   Ode_Jacobian        = @arenstorfOrbit_Jacobian;
%   Ode_YTLM            = eye(4);
%   Time_Interval       = [ 0 17.0652166 ];
%   Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];
%
% Now that we have our model loaded in our workspace, we can perform a
% tangent linear explicit Runge-Kutta integration using MATLODE's prebuilt default
% settings. We note that a Jacobian and Y_TLM are required for sensitivity analysis.
%
%   Options  = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Y_TLM',Ode_YTLM);
%   [ T, Y, Sens ] = MATLODE_SDIRK_TLM_Integrator(Ode_Function,Time_Interval,Y0,Options);
%
% Printing out our results, we can analyze our model state at our final
% time.
%
%   disp('solution and tangent linear sensitivty at Time_Interval(2)');
%   disp(Y);
%   disp(Sens);
%
% For addition examples, see Help -> Supplemental Software -> Examples ->
% Sensitivity Analysis -> MATLODE_Example_SDIRK_TLM_Integrator.
%
%% Contact Information
%%
% Dr. Adrian Sandu                 | Phone: (540) 231-2193 | Email: sandu@cs.vt.edu
%%
% Tony D'Augustine                 | Phone: (540) 231-6186 | Email: adaug13@vt.edu 
%%
% Computational Science Laboratory | Phone: (540) 231-6186
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
% [2] Hong Zhang, Adrian Sandu. FATODE: a library for forward, adjoint and 
%     tangent linear integration of ODEs, SIAM Journal on Scientific 
%     Computing, 36(5), C504–C523, 2014
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ Tout_TLM, Yout_TLM, Y_TLM, Stats_TLM ] = MATLODE_SDIRK_TLM_Integrator( OdeFunction, Tspan, Y0, OPTIONS_U )
    % Display input/output parameters
    if ( nargout == 0 && nargin == 0 )
        fprintf('Syntax: \n');
        fprintf( '                          MATLODE_SDIRK_TLM_Integrator \n' );
        fprintf( '         [ T, Y, Sens ] = MATLODE_SDIRK_TLM_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
        fprintf( '  [ T, Y, Sens, Stats ] = MATLODE_SDIRK_TLM_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
        fprintf( '\n' );
        return;
    end 

    % Force initial value matrix to be N X 1.
    if ( size(Y0,2) == 1 )
        % DO NOTHING
    else
        Y0 = transpose(Y0);
    end     
    
    % Configure Options
    [ OPTIONS, Coefficient ] = OPTIONS_Configuration(OPTIONS_U,'SDIRK','TLM',Y0, Tspan );
    
    % Check input dimensions
    OPTIONS = Input_Dimension(Tspan(1), Y0, OdeFunction, OPTIONS);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call Tangent Linear Method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Accounts for tspan n-dimensional array to force output at
    %   specific times. Accumlates statistics.
    tspanMaxSize = max(size(Tspan));
    Yout_TLM_interval = transpose(Y0);
    Yout_TLM = transpose(Y0);
    Tout_TLM = Tspan(1);
    ISTATUS_TLM = ISTATUS_Struct('default');
    for interval=1:tspanMaxSize-1
        tic;
        [ Tout_TLM_interval, Yout_TLM_interval, Y_TLM_Interval, ISTATUS_TLM_interval, RSTATUS_TLM, Ierr ] = ...
            SDIRK_TLM_Integrator( OdeFunction,[Tspan(interval), Tspan(interval+1)], Yout_TLM_interval(end,:), OPTIONS, Coefficient);
        OPTIONS.Y_TLM = Y_TLM_Interval;
        elapsedTime(interval) = toc;
        ISTATUS_TLM = ISTATUS_Add(ISTATUS_TLM,ISTATUS_TLM_interval);
        Tout_TLM = [Tout_TLM; Tout_TLM_interval];
        Yout_TLM = [Yout_TLM; Yout_TLM_interval];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Fold statistics
    RSTATUS_TLM.Etime = sum(elapsedTime);
    Stats_TLM.ISTATUS = ISTATUS_TLM;
    Stats_TLM.RSTATUS = RSTATUS_TLM;    
    
    % Display Statistics
    if ( OPTIONS.displayStats == true )
        fprintf( '\n SDIRK_TLM Statistics' );
        PrintISTATUS( ISTATUS_TLM);
        PrintRSTATUS( RSTATUS_TLM );
    end
    
    % Outputs
    Y_TLM = OPTIONS.Y_TLM;
        
return;

%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt,edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
%
%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>