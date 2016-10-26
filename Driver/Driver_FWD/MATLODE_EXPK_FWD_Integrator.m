%% MATLODE_EXPK_FWD_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%                   MATLODE_EXPK_FWD_Integrator
%          [T, Y] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)
%   [T, Y, Stats] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)
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
% Driver file to solve the system y' = F(t,y) using a Exponential-Krylov (EXPK)
% method.
%
% |MATLODE_EXPK_FWD_Integrator| displays the available methods
% associated with the exponential forward integrator.
%
% |[T, Y] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options)|
% computes the ODE solution with respect to the user supplied options
% configuration.
%
% |[T, Y, Stats] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0,
% Options)| computes the ODE solution with respect to the user supplied
% options configuration and returns the computation statistics.
%
%% Example
% For the following examples we will use Arenstorf Orbit as a toy problem 
% to illustrate |MATLODE_EXP_FWD_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load our input parameters into our workspace. 
%
%   Ode_Function        = @arenstorfOrbit_Function;
%   Ode_Jacobian        = @arenstorfOrbit_Jacobian;
%   Time_Interval       = [ 0 17.0652166 ];
%   Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];
%
% Now that we have our model loaded in our workspace, we can perform a
% forward exponential integration using MATLODE's prebuilt default
% settings. We note that a Jacobian is required for an expoenential-krylov method.
%
%   Options  = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'AbsTol',1d-7,'RelTol',1d-7);
%   [ T, Y ] = MATLODE_EXPK_FWD_Integrator(Ode_Function,Time_Interval,Y0,Options);
%
% Printing out our results, we can analyze our model state at our final
% time.
%
%   disp('solution at Time_Interval(2)');
%   disp(Y);
%
% For addition examples, see Help -> Supplemental Software -> Examples ->
% Forward Integration -> MATLODE_Example_EXPK_FWD_Integrator.
%
%% Contact Information
%%
% Dr. Adrian Sandu                 | Phone: (540) 231-2193 | Email: sandu@cs.vt.edu
%%
% Paul Tranquilli                  | Phone: (540) 231-6186 | Email: ptranq@vt.edu 
%%
% Computational Science Laboratory | Phone: (540) 231-6186
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE
%
% [2] P. Tranquilli, A. Sandu, 'Exponential-Krylov methods for ordinary differential equations,' Journal of Computational Physics, 278, 31-46
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function  [ Tout_FWD, Yout_FWD, Stats_FWD ] = MATLODE_EXPK_FWD_Integrator( OdeFunction, Tspan, Y0, OPTIONS_U )
    % Display input/output parameters
    if ( nargout == 0 && nargin == 0 )
        fprintf('Syntax: \n');
        fprintf( '                     MATLODE_EXPK_FWD_Integrator \n' );
        fprintf( '          [ T, Y ] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
        fprintf( '   [ T, Y, Stats ] = MATLODE_EXPK_FWD_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
        fprintf( '\n' );
        return;
    end
    
    % Force initial value matrix to be 1 X N.
    if ( size(Y0,2) == 1 )
        % DO NOTHING
    else
        Y0 = transpose(Y0);
    end      
    
    % Initialize Adjoint Flags
    adjStackFlag = false;
    adjQuadFlag  = false;    
    
    % Configure Options
    [ OPTIONS, Coefficient ] = OPTIONS_Configuration(OPTIONS_U,'EXPK','FWD',Y0, Tspan);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call Forwrad Method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Accounts for tspan n-dimensional array to force output at
    %   specific times. Accumlates statistics.
    tspanMaxSize = max(size(Tspan));
    Yout_FWD_interval = Y0';
    Yout_FWD = Y0';
    Tout_FWD = Tspan(1);
    ISTATUS_FWD = ISTATUS_Struct('default');
    for interval=1:tspanMaxSize-1
        tic;
        [ Tout_FWD_interval, Yout_FWD_interval, ISTATUS_FWD_interval, RSTATUS_FWD, Ierr ] = ...
            EXPK_FWD_Integrator( OdeFunction,[Tspan(interval), Tspan(interval+1)], Yout_FWD_interval(end,:), OPTIONS, Coefficient, adjStackFlag, adjQuadFlag );
        elapsedTime(interval) = toc;
        ISTATUS_FWD = ISTATUS_Add(ISTATUS_FWD,ISTATUS_FWD_interval);
        Tout_FWD = [Tout_FWD; Tout_FWD_interval];
        Yout_FWD = [Yout_FWD; Yout_FWD_interval];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fold statistics
    RSTATUS_FWD.Etime = sum(elapsedTime);
    Stats_FWD.ISTATUS = ISTATUS_FWD;
    Stats_FWD.RSTATUS = RSTATUS_FWD;     
    
    % Display Statistics
    if ( OPTIONS.displayStats == true )
        fprintf( '\nForward Statistics' );
        PrintISTATUS( ISTATUS_FWD );
        PrintRSTATUS( RSTATUS_FWD );
    end

return
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
%       <td>Paul Tranquilli</td>
%       <td>ptranq@vt.edu</td>
%       <td>Implemented ROK and all supporting components</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13#vt.edu</td>
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