%% MATLODE_ERK_FWD_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%                   MATLODE_ERK_FWD_Integrator
%          [T, Y] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0)
%   [T, Y, Stats] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0)
%          [T, Y] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)
%   [T, Y, Stats] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)
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
% |Stats|: integrator statistics
%
%% Description
% Driver file to solve the system y' = F(t,y) using a Explicit Runge-Kutta (ERK)
% method.
%
% |MATLODE_ERK_FWD_Integrator| displays the available methods
% associated with the explicit Runge Kutta forward integrator.
%
% |[T, Y] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0)| computes the ODE
% solution at the final time using the default parameters.
%
% |[T, Y, Stats] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0)|
% computes the ODE solution at the final time using the default
% parameters and returns computational statistics.
%
% |[T, Y] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)|
% computes the ODE solution with respect to the user supplied options
% configuration.
%
% |[T, Y, Stats] = MATLODE_ERK_FWD_Integrator(Ode_Function, Time_Interval, Y0,
% Options)| computes the ODE solution with respect to the user supplied
% options configuration and returns the computation statistics.
%
%% Example
% For the following examples we will use Arenstorf Orbit as a toy problem 
% to illustrate |MATLODE_ERK_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load our input parameters into our workspace. 
%
%   Ode_Function        = @arenstorfOrbit_Function;
%   Time_Interval       = [ 0 17.0652166 ];
%   Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];
%
% Now that we have our model loaded in our workspace, we can perform a
% forward explicit Runge-Kutta integration using MATLODE's prebuilt default
% settings.
%
%   [ T, Y ] = MATLODE_ERK_FWD_Integrator(Ode_Function,Time_Interval,Y0);
%
% Printing out our results, we can analyze our model state at our final
% time.
%
%   disp('solution at Time_Interval(2)');
%   disp(Y);
%
% For addition examples, see Help -> Supplemental Software -> Examples ->
% Forward Integration -> MATLODE_Example_ERK_FWD_Integrator.
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
%     Computing, 36(5), C504-C523, 2014.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ Tout_FWD, Yout_FWD, FWD_Stats ] = ...
    MATLODE_EXIM_MRGARK_FWD_Integrator( OdeFunctionOne, OdeFunctionTwo, Tspan, Y0, varargin )

    % Display input/output parameters
    if ( nargout == 0 && nargin == 0 )
        fprintf('Syntax: \n');
        fprintf( '                     MATLODE_MRGARK_ERK_FWD_Integrator \n' );
        fprintf( '          [ T, Y ] = MATLODE_MRGARK_ERK_FWD_Integrator(odefunone, odefuntwo, Tspan, y0) \n' );
        fprintf( '   [ T, Y, Stats ] = MATLODE_MRGARK_ERK_FWD_Integrator(odefunone, odefuntwo, Tspan, y0) \n' );
        fprintf( '          [ T, Y ] = MATLODE_MRGARK_ERK_FWD_Integrator(odefunone, odefuntwo, Tspan, y0, options) \n' );
        fprintf( '   [ T, Y, Stats ] = MATLODE_MRGARK_ERK_FWD_Integrator(odefunone, odefuntwo, Tspan, y0, options) \n' );
        fprintf( '\n' );
        return;
    end

    % Check inputs
    if ( nargin < 4 )
        error('Not enough input arguments.');
    elseif ( nargin == 4 )
        OPTIONS_U = MATLODE_OPTIONS();
    elseif ( nargin == 5 )
        OPTIONS_U = varargin{1};
    else
        error('Too many input arguments.');
    end
    
    % Initialize Adjoint Flags
    adjQuadFlag = false;
        
    % Configure Forward Options
    [ OPTIONS, Coefficient ] = OPTIONS_Configuration(OPTIONS_U, 'EXIMMRGARK', 'FWD', Y0, Tspan );
    
    % Check input dimensions
    OPTIONS = Input_Dimension(Tspan, Y0, {OdeFunctionOne, OdeFunctionTwo}, OPTIONS);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call Forwrad Method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Accounts for Tspan n-dimensional array to force output at
    %   specific times. Accumlates statistics.
    tspanMaxSize = max(size(Tspan));
    Yout_FWD = transpose(Y0);
    Tout_FWD = Tspan(1);
    ISTATUS_FWD = ISTATUS_Struct('default');

    for interval=1:tspanMaxSize-1
        tic;
        [ Tout_FWD_interval, Yout_FWD_interval, ISTATUS_FWD_interval, RSTATUS_FWD, Ierr ] = ...
            EXIM_MRGARK_FWD_Integrator( OdeFunctionOne, OdeFunctionTwo, [Tspan(interval), Tspan(interval+1)], Y0, OPTIONS, ...
            Coefficient, adjQuadFlag );
        elapsedTime(interval) = toc;
        ISTATUS_FWD = ISTATUS_Add(ISTATUS_FWD,ISTATUS_FWD_interval);
         Tout_FWD = [Tout_FWD; transpose(Tout_FWD_interval)];
         Yout_FWD = [Yout_FWD;Yout_FWD_interval];
        Y0=transpose(Yout_FWD(end,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fold statistics
    RSTATUS_FWD.Etime = sum(elapsedTime);
    FWD_Stats.ISTATUS = ISTATUS_FWD;
    FWD_Stats.RSTATUS = RSTATUS_FWD;
    
% %    Display Statistics
%     if ( OPTIONS.displayStats == true )
%        fprintf( '\nForward Statistics' );
%        PrintISTATUS( ISTATUS_FWD );
%        PrintRSTATUS( RSTATUS_FWD );
%     end
    
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