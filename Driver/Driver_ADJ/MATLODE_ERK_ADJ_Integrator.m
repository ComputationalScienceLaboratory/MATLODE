%% MATLODE_ERK_ADJ_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%                                     MATLODE_ERK_ADJ_Integrator
%                    [ T, Y, Sens ] = MATLODE_ERK_ADJ_Integrator(Ode_Function, Time_Interval, Y0, Options)
%   [ T, Y, Sens, Quad, Mu, Stats ] = MATLODE_ERK_ADJ_Integrator(Ode_Function, Time_Interval, Y0, Options)
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
% |Mu|: Mu term
%
% |Stats|: integrator statistics
%
%% Description
% Driver file to solve the system y' = F(t,y) and adjoint sensitivity
% using an Explicit Runge-Kutta (ERK) method.
%
% |MATLODE_ERK_ADJ_Integrator| displays the available methods
% associated with the exponential forward integrator.
%
% |[T, Y, Sens] = MATLODE_ERK_ADJ_Integrator(Ode_Function, Time_Interval, Y0,
% Options)| computes the ODE solution with respect to the user supplied
% options configuration and adjoint sensitivity.
%
% |[T, Y, Sens, Quad, Mu, Stats]| = MATLODE_ERK_ADJ_Integrator(Ode_Function,
% Time_Interval, Y0, Options)| computes the ODE solution with respect to
% the user supplied options configuration, sensitivity, quadrature, mu and
% statisitics.
%
%% Example
% For the following examples we will use Arenstorf Orbit as a toy problem
% to illustrate |MATLODE_ERK_ADJ_Integrator| functionalities and features.
% To initially setup Arenstorf Orbit, execute the MATLAB commands below to
% load our input parameters into our workspace.
%
%   Ode_Function        = @arenstorfOrbit_Function;
%   Ode_Jacobian        = @arenstorfOrbit_Jacobian;
%   Ode_Lambda          = eye(4);
%   Time_Interval       = [ 0 17.0652166 ];
%   Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];
%
% Now that we have our model loaded in our workspace, we can perform a
% adjoint explicit Runge-Kutta integration using MATLODE's prebuilt default
% settings. We note that a Jacobian and Lambda are required for sensitivity analysis.
%
%   Options  = MATLODE_OPTIONS('Jacobian',Ode_Jacobian,'Lambda',Ode_Lambda);
%   [ T, Y, Sens ] = MATLODE_ERK_ADJ_Integrator(Ode_Function,Time_Interval,Y0,Options);
%
% Printing out our results, we can analyze our model state at our final
% time.
%
%   disp('solution and adjoint sensitivty at Time_Interval(2)');
%   disp(Y);
%   disp(Sens);
%
% For addition examples, see Help -> Supplemental Software -> Examples ->
% Sensitivity Analysis -> MATLODE_Example_ERK_ADJ_Integrator.
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
function [ Tout_FWD, Yout_FWD, Lambda, Quadrature, Mu, Stats ] = ...
    MATLODE_ERK_ADJ_Integrator( OdeFunction, Tspan, Y0, OPTIONS_U )
% Display input/output parameters
if ( nargout == 0 && nargin == 0 )
    fprintf('Syntax: \n');
    fprintf( '                                  MATLODE_ERK_ADJ_Integrator \n' );
    fprintf( '                 [ T, Y, Sens ] = MATLODE_ERK_ADJ_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
    fprintf( '[ T, Y, Sens, Quad, Mu, Stats ] = MATLODE_ERK_ADJ_Integrator(Ode_Function, Time_Interval, Y0, Options) \n' );
    fprintf( '\n' );
    return;
end

% Initialize
stack_ptr    = 0;
% Adjoint Flags
adjStackFlag = true;
adjQuadFlag  = false;
adjMuFlag    = false;

% Initialize ISTATUS and RSTATUS
ISTATUS_TOTAL = ISTATUS_Struct('default');

[ OPTIONS, Coefficient ] = OPTIONS_Configuration(OPTIONS_U, 'ERK', 'ADJ', Y0, Tspan );
OPTIONS                  = Input_Dimension(Tspan(1), Y0, OdeFunction, OPTIONS);

% Adjoint Flags
if ( ~isempty(OPTIONS.Jacp) && ~isempty(OPTIONS.Mu) )
    if ( ~isempty(OPTIONS.QFun) && ~isempty(OPTIONS.DRDP) && ...
            ~isempty(OPTIONS.DRDY) && ~isempty(OPTIONS.Quadrature) )
        
        adjQuadFlag  = true;
        adjMuFlag    = true;
    else
        
        adjQuadFlag  = false;
        adjMuFlag    = true;
    end
else
    
    if ( ~isempty(OPTIONS.QFun) && ~isempty(OPTIONS.DRDY) && ~isempty(OPTIONS.Quadrature) )
        
        adjQuadFlag  = true;
        adjMuFlag    = false;
    else
        
        adjQuadFlag  = false;
        adjMuFlag    = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call Forward and Adjoint Core Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Accounts for tspan n-dimensional array to force output at
%   specific times. Accumlates statistics.

tspanMaxSize = max(size(Tspan));
Yout_FWD = transpose(Y0);
Tout_FWD = Tspan(1);
FWD_ISTATUS = ISTATUS_Struct('default');
for interval=1:tspanMaxSize-1
    tic;
    [ FWD_Tout_Interval, FWD_Yout_Interval, FWD_ISTATUS_interval, FWD_RSTATUS, FWD_Ierr, stack_ptr, Quadrature ] = ...
        ERK_FWD_Integrator( OdeFunction,[Tspan(interval), Tspan(interval+1)], Y0, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr );
    elapsedTime_FWD(interval) = toc;
    FWD_ISTATUS = ISTATUS_Add(FWD_ISTATUS,FWD_ISTATUS_interval);
    Tout_FWD = [Tout_FWD; transpose(FWD_Tout_Interval)];
    Yout_FWD = [Yout_FWD; FWD_Yout_Interval];
    Y0=transpose(Yout_FWD(end,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call ERK ADJ Core Method ;
Tf = Tout_FWD(end);
Yf = Yout_FWD(end,:);
Lambda_Tf = OPTIONS.Lambda(Tf,Yf);
if adjMuFlag
    Mu_Tf = OPTIONS.Mu(Tf,Yf);
else
    Mu_Tf = NaN;
end

tic;
[ ADJ_Tout, ADJ_Yout, Lambda, Mu, ADJ_ISTATUS, ADJ_RSTATUS, ADJ_Ierr ] = ...
    ERK_ADJ_DiscreteIntegrator(Lambda_Tf,Mu_Tf,OPTIONS, Coefficient, stack_ptr, adjQuadFlag,adjMuFlag );
elapsedTime_ADJ = toc;

% Fold statistics
FWD_RSTATUS.Etime = sum(elapsedTime_FWD);
ADJ_RSTATUS.Etime = elapsedTime_ADJ;
Stats.ISTATUS_FWD = FWD_ISTATUS;
Stats.RSTATUS_FWD = FWD_RSTATUS;
Stats.ISTATUS_ADJ = ADJ_ISTATUS;
Stats.RSTATUS_ADJ = ADJ_RSTATUS;

% Display Forward Statistics
if ( OPTIONS.displayStats == true )
    fprintf( '\nForward Statistics' );
    PrintISTATUS( FWD_ISTATUS );
    PrintRSTATUS( FWD_RSTATUS );
end

% Display Adjoint Statistics
if ( OPTIONS.displayStats == true )
    fprintf( '\n Adjoint Statistics' );
    PrintISTATUS( ADJ_ISTATUS );
    PrintRSTATUS( ADJ_RSTATUS );
end

ISTATUS_TOTAL.Nfun = FWD_ISTATUS.Nfun + ADJ_ISTATUS.Nfun;
ISTATUS_TOTAL.Njac = FWD_ISTATUS.Njac + ADJ_ISTATUS.Njac;
ISTATUS_TOTAL.Nstp = FWD_ISTATUS.Nstp + ADJ_ISTATUS.Nstp;
ISTATUS_TOTAL.Nacc = FWD_ISTATUS.Nacc + ADJ_ISTATUS.Nacc;
ISTATUS_TOTAL.Nrej = FWD_ISTATUS.Nrej + ADJ_ISTATUS.Nrej;
ISTATUS_TOTAL.Ndec = FWD_ISTATUS.Ndec + ADJ_ISTATUS.Ndec;
ISTATUS_TOTAL.Nsol = FWD_ISTATUS.Nsol + ADJ_ISTATUS.Nsol;
ISTATUS_TOTAL.Nsng = FWD_ISTATUS.Nsng + ADJ_ISTATUS.Nsng;
ISTATUS_TOTAL.Nchk = FWD_ISTATUS.Nchk + ADJ_ISTATUS.Nchk;

% Display Total Statistics
if ( OPTIONS.displayStats == true )
    fprintf( '\n Total Statistics' );
    PrintISTATUS( ISTATUS_TOTAL );
end

% Outputs
if ( ~exist('Quadrature','var') )
    Quadrature = [];
end
if ( ~exist('Mu','var')  )
    Mu = [];
end


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