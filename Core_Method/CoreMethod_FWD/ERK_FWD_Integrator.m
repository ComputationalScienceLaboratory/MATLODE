%% ERK_FWD_Integrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
% <html> Up: <a href="../../../Library/html/Library.html">Library</a> </html>
%
%% Syntax
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ERK_FWD_Integrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )
%
%% Input Parameters
% |OdeFunction|: ODE function function handle
%
% |Tspan|: Time interval
%
% |OPTIONS|: Option struct
%
% |Coefficients|: Constant coefficients associated with method
%
% |adjStackFlag|: Adjoint snapshot stack flag
%
% |adjQuadFlag|: Adjoint quadrature stack flag
%
% |stack_ptr|: pointer for global snapshot stack
%
%% Output Parameters
% |Tout|: Time vector
%
% |Yout|: State vector
%
% |ISTATUS|: Integer statistics
%
% |RSTATUS|: Real statistics
%
% |Ierr|: Error flag
%
% |stack_ptr|: Pointer for global snapshot stack
%
% |quadrature|: Adjoint quadrature
%
%% Description
% Explicit Runge-Kutta forward core method.
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
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ERK_FWD_Integrator( OdeFunction,...
    Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global rkA rkB rkC rkE

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get Problem Size
NVAR = OPTIONS.NVAR;

Tinitial = Tspan(1);
Tfinal = Tspan(2);

Roundoff = eps/2;

% Allocate Memory
[ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );

TYindex = 1;



Yout(:,TYindex) = Y;
Tout(TYindex,1) = Tinitial;
if ( adjQuadFlag == true ) % initial value of the quadrature
    quadrature = OPTIONS.Quadrature(Tinitial,Y);
else
    quadrature = NaN;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ISTATUS = ISTATUS_Struct('default');
RSTATUS = RSTATUS_Struct('default');

ISTATUS.Nchk = 1;




T = Tinitial;
Tdirection = sign( Tfinal-Tinitial );
H = max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) );

if ( abs(H) <= 10.0*Roundoff)
    H = 1.0e-6;
end

H = min( abs(H), OPTIONS.Hmax );
sign_Tdirection = sign(Tdirection );
H = H*sign_Tdirection;
Reject = false;
FirstStep = true;

% Determine scaling factor for integration
SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while ( (Tfinal-T) * Tdirection - Roundoff >= 0.0 )
    if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
        error('Number of steps exceeds maximum bound.');
    end
    
    % case 1: determines if we are close enough to tfinal due to
    %         roundoff error
    % case 2: cannot distinguish between tnext and tcurrent
    if ( abs( Tfinal-T) <= 10*Roundoff*abs(Tfinal))
        T = Tfinal;
        break;
    elseif ( (T+0.1*H == T) || (abs(H) <= Roundoff) )
        error('Step size too small; T + 10*H = T or H < eps/2.');
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Stages
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    K = zeros(NVAR,Coefficient.NStage);
    Z = zeros(NVAR,Coefficient.NStage);
    for istage = 1:Coefficient.NStage
        if ( istage > 1 )
            for j = 1:istage-1
                Z(:,istage) = Z(:,istage) + H*rkA(istage, j)*K(:,j);
            end
        end
        Z(:,istage) = Z(:,istage) + Y;
        K(:,istage) = OdeFunction(T+rkC(istage)*H,Z(:,istage));
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Error estimation
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS.Nstp = ISTATUS.Nstp + 1;
    TMP = zeros(NVAR,1);
    for i = 1:Coefficient.NStage
        if (rkE(i) ~= 0.0 )
            TMP = TMP + H*rkE(i)*K(:,i);
        end
    end
    
    % Perform error norm
    Err = sum((TMP.*SCAL).^2,1);
    Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );
    
    % Computation of new step size Hnew
    Fac = OPTIONS.FacSafe*(Err)^(-1.0/Coefficient.ELO);
    Fac = max( OPTIONS.FacMin, min(OPTIONS.FacMax, Fac) );
    Hnew = H*Fac;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Accept/Reject step
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( Err < 1.0 )
        FirstStep = false;
        ISTATUS.Nacc = ISTATUS.Nacc + 1;
        
        % Checkpoints for Adjoint Calculations
        if ( adjStackFlag == true )
            stack_ptr = ERK_Push( NVAR, Coefficient.NStage, T, H, Y, Z, OPTIONS.Max_no_steps, stack_ptr );
        end
        
        % Update the results for the quadrature term
        if ( adjQuadFlag == true )
            for i=1:Coefficient.NStage
                TMP =  Z(:,i);
                R = OPTIONS.QFun(T+rkC(i)*H,TMP);
                quadrature = H*rkB(i)*R + quadrature;
            end
        end
        
        % Update Memory Allocation
        if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
            [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
            ISTATUS.Nchk = ISTATUS.Nchk + 1;
        end
        
        % Update time
        T = T + H;
        
        % Update solution
        for i = 1:Coefficient.NStage
            if ( rkB(i) ~= 0.0 )
                Y = Y + H*rkB(i)*K(:,i);
            end
        end
        
        % Store checkpoint values
        if ( OPTIONS.storeCheckpoint == true )
            Tout(TYindex,1) = T;
            Yout(:,TYindex) = Y;
            TYindex = TYindex + 1;
        end
        
        % Update scaling coefficients
        SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
        
        % Next time step
        Hnew = Tdirection*min( abs(Hnew), OPTIONS.Hmax);
        
        % Last T and H
        RSTATUS.Ntexit = T;
        RSTATUS.Nhexit = H;
        RSTATUS.Nhnew = Hnew;
        
        % No step increase after a rejection
        if ( Reject )
            Hnew = Tdirection*min( abs(Hnew), abs(H) );
        end
        Reject = false;
        
        % Update step size
        if ( (T+Hnew-Tfinal)*Tdirection > 0.0)
            H = Tfinal - T;
        else
            H = max(OPTIONS.Hmin, min(OPTIONS.Hmax, Hnew));
        end
        
        % for debugging
        if ( OPTIONS.displaySteps == true )
            str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
            disp(str);
        end
        
    else % Step is rejected
        % Update step size
        if ( FirstStep || Reject )
            H = OPTIONS.FacRej*H;
        else
            H = max(OPTIONS.Hmin, min(OPTIONS.Hmax, Hnew));
        end
        Reject = true;
        
        % Update number of rejections
        if ( ISTATUS.Nacc >= 1.0 )
            ISTATUS.Nrej = ISTATUS.Nrej + 1;
        end
        
        %for debugging
        if ( OPTIONS.displaySteps == true )
            str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H) ];
            disp(str);
        end
        
    end % END accept
end % Time loop

% Successful return
Ierr = 1;

% Deallocate Memory
if ( OPTIONS.storeCheckpoint == true )
    %[ Tout, Yout ] = matlOde_deallocateMemory(Tout,Yout,ISTATUS.Nacc);
    Tout(ISTATUS.Nacc+1:end,:) = [];
    Yout(:,ISTATUS.Nacc+1:end) = [];
else
    Tout = T;
    Yout = Y;
end
Tout = transpose(Tout);
Yout = transpose(Yout);

return;



% (START) SUBROUTINE: Explcit Runge Kutta Push
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Saves the next trajectory snapshot for discrete adjoints
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function stack_ptr = ERK_Push( NVAR, rkS, T, H, Y, Z, Max_no_steps, stack_ptr )
global chk_H chk_T chk_Y chk_Z

stack_ptr = stack_ptr + 1;
if ( stack_ptr > Max_no_steps )
    error( 'Push failed: buffer overflow' );
end
chk_H( stack_ptr ) = H;
chk_T( stack_ptr ) = T;
chk_Y( 1:NVAR, stack_ptr ) = Y(1:NVAR);
chk_Z( 1:NVAR, 1:rkS, stack_ptr ) = Z(1:NVAR, 1:rkS);
%chk_Z( 1:rkS, 1:NVAR, stack_ptr ) = Z;

return; % (END) SUBROUTINE: Explicit Runge Kutta Push

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
