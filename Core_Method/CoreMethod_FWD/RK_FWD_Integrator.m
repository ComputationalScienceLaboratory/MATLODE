%% RK_FWD_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = RK_FWD_Integrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )
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
% Implicit Runge-Kutta forward core method.
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
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = RK_FWD_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkTheta
    global rkBgam
    global rkGamma
    global rkB rkC rkD rkF
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end 

    % Get Problem Size
    NVAR = max(size(Y));    

    Tinitial = Tspan(1);
    Tfinal = Tspan(end); 
    
    T = Tinitial;

    Yout = zeros(NVAR,OPTIONS.Max_no_steps);
    Tout = zeros(OPTIONS.Max_no_steps,1);
    TYindex = 1;
    
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;     
    
    Roundoff = 1.11022302462515654E-016;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
    
    ISING = false;
    Transp = false;
        
    CONT = zeros(NVAR,1);
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial Setting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    quadrature = OPTIONS.Quadrature();

    Tdirection = sign( Tfinal-Tinitial );
    H = min( max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) ), OPTIONS.Hmax );

    if ( abs(H) <= 10.0*Roundoff ) 
        H = 1.0d-6;
    end

    H = sign(Tdirection)*H;
    Hold = H;
    Reject    = false;
    FirstStep = true;
    SkipJac   = false;
    SkipLU    = false;

    if ( ((Tinitial + H*1.0001 - Tfinal)*Tdirection) >= 0.0 )
        H = Tfinal - Tinitial;
    end

    Nconsecutive = 0;
    
    % Determine scaling factor for integration
    SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( (Tfinal-T)*Tdirection - Roundoff > 0.0 ) % Tloop
        
        if ( ~SkipLU )
            % Compute the jacobian matrix
            if ( ~SkipJac )
                JAC = OPTIONS.Jacobian(T,Y);
                ISTATUS.Njac = ISTATUS.Njac + 1;
            end

            % Compute the matrices E1 and E2 and their decompositions
            [ e1, e2, ISING, ISTATUS ] = RK_Decomp( NVAR, H, JAC, ISING, ISTATUS );
            if ( ISING ~= 0 )
                ISTATUS.Nsng = ISTATUS.Nsng + 1;
                Nconsecutive = Nconsecutive + 1;
                if ( Nconsecutive >= 5 )
                    error('Non-convergence of Newton iterations.');
                end
                H = H*0.5;
                Reject  = true; 
                SkipJac = true; 
                SkipLU  = false;
                % CYCLE Tloop <---- NEED TO LOOK UP WHAT CYCLE DOES IN F90
            else
                Nconsecutive = 0;
            end
        end

        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
            str = ['Max number of time steps is ', num2str(OPTIONS.Max_no_steps)];
            disp(str);
            error('Number of steps exceeds maximum bound.');
        end

        if ( 0.1*abs(H) <= abs(T)*Roundoff )
            error('Step size too small: T + 10*H = T or H < Roundoff');
        end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loop for simplified Newton iterations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Starting values for Newton iteration
        if ( FirstStep || ( ~OPTIONS.StartNewton ) )
            Z1 = zeros(NVAR,1);
            Z2 = zeros(NVAR,1);
            Z3 = zeros(NVAR,1);
        else
            % Evaluate quadratic polynomial
            [ Z1, Z2, Z3, CONT ] = RK_Interpolate( 'eval', NVAR, H, Hold, Z1, Z2, Z3, CONT );
        end

        % Initializations for Newton iteration
        NewtonDone = false;
        Fac = 0.5; % Step reduction if too many iterations

        for NewtonIter=1:OPTIONS.NewtonMaxIt % NewtonLoop
                        
            % Prepare the right-hand side
            [ DZ1, DZ2, DZ3, F0, L3A_flag ] = RK_PrepareRHS( T, H, Y, Z1, Z2, Z3, OdeFunction, Coefficient );
             
            % Solve for linear systems
            [ DZ1, DZ2, DZ3, ISTATUS ] = RK_Solve( NVAR, H, DZ1, DZ2, DZ3, ...
                e1, e2, ISING, ISTATUS );

            NewtonIncrement = sqrt( ( errorNorm( NVAR, DZ1, SCAL )^2 ...
                + errorNorm( NVAR, DZ2, SCAL )^2 ...
                + errorNorm( NVAR, DZ3, SCAL )^2 )/3.0 );
             
            if ( NewtonIter == 1 ) 
                Theta = abs(OPTIONS.ThetaMin);
                NewtonRate = 2.0;
            else
                Theta = NewtonIncrement/NewtonIncrementOld;
                if ( Theta < 0.99 ) 
                    NewtonRate = Theta/( 1.0-Theta );
                else % Non-convergence of Newton: Theta too large
                    break; % NewtonLoop
                end
                if ( NewtonIter < OPTIONS.NewtonMaxIt ) 
                    % Predict error at the end of Newton process
                    NewtonPredictedErr = NewtonIncrement*Theta^( OPTIONS.NewtonMaxIt ...
                        - NewtonIter )/( 1.0-Theta );
                    if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                        % Non-convergence of Newton: predicted error too large
                        Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS.NewtonTol );
                        Fac = 0.8*Qnewton^( -1.0/( 1+OPTIONS.NewtonMaxIt - NewtonIter ));
                        break; % NewtonLoop
                    end
                end
            end

            NewtonIncrementOld = max( NewtonIncrement, Roundoff );
            % Update Solution
            Z1 = Z1 - DZ1;
            Z2 = Z2 - DZ2;
            Z3 = Z3 - DZ3;
            
            % Check error in Newton iterations
            NewtonDone = ( NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
                        
            if ( NewtonDone )
                NewIt = NewtonIter;
                break; % NewtonLoop
            end
            if ( NewtonIter == OPTIONS.NewtonMaxIt ) 
                disp( 'Slow or no convergence in Newton Iteration: Max no. of Newton iterations reached' );
            end

        end % NewtonLoop

        if ( ~NewtonDone )
            H       = Fac*H;
            Reject  = true;
            SkipJac = true;
            SkipLU  = false;
            continue; % Tloop
        end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SDIRK Stage
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( OPTIONS.SdirkError )
            % Starting values for Newton iterations
            Z4 = Z3;
            
            % Prepare the loop-independent part of the right-hand side
            % G = H*rkBgam(0)*F0 + rkTheta(1)*Z1 + rkTheta(2)*Z2 + rkTheta(3)*Z3
            G = zeros(NVAR,1);
            % L3A = 5
            if ( Coefficient.Method ~= 5 )
                if ( L3A_flag == 0 )
                    F0 = OdeFunction(T,Y);
                    ISTATUS.Nfun = ISTATUS.Nfun + 1;
                end
                G = G + rkBgam(1)*H*F0;
            end
            G = G + rkTheta(1)*Z1;
            G = G + rkTheta(2)*Z2;
            G = G + rkTheta(3)*Z3;
            
            % Initializations for Newton iteration
            NewtonDone = false;
            Fac = 0.5; % Step reduction factor if too many iterations

            for NewtonIter = 1:OPTIONS.NewtonMaxIt % SDNewtonLoop
                % Prepare the loop-dependent part of the right-hand side
                TMP = Y + Z4;
                DZ4 = OdeFunction( T+H,TMP );
                ISTATUS.Nfun = ISTATUS.Nfun + 1;
                % DZ4(1:NVAR) = { G(1:NVAR) - Z4(!:NVAR)*(rkGamma/H) + DZ4(1:NVAR)
                DZ4 = DZ4 - (rkGamma*Z4)/H;
                DZ4 = DZ4 + (rkGamma*G)/H;

                % Solve the linear system              
                DZ4 = e1\DZ4;
                
                ISTATUS.Nsol = ISTATUS.Nsol + 1;

                % Check convergence of Newton iterations
                NewtonIncrement = errorNorm( NVAR, DZ4, SCAL );
                if ( NewtonIter == 1 ) 
                    ThetaSD = abs( OPTIONS.ThetaMin );
                    NewtonRate = 2.0;
                else
                    ThetaSD = NewtonIncrement/NewtonIncrementOld;
                    if ( ThetaSD < 0.99 ) 
                        NewtonRate = ThetaSD / (1.0 - ThetaSD);
                        % Predict error at the end of Newton process
                        NewtonPredictedErr = NewtonIncrement*ThetaSD^( OPTIONS.NewtonMaxIt ...
                            - NewtonIter )/( 1.0 - ThetaSD );
                        if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                            % Non-convergence of Newton: predicted error too large
                            % str = ['Error too large: ', NewtonPredictedErr ];
                            % disp(str);
                            Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS.NewtonTol );
                            Fac = 0.8*Qnewton^( -1.0/( 1.0 + OPTIONS.NewtonMaxIt - NewtonIter ) );
                            break;
                            % EXIT SDNewtonLoopl
                        end
                    else % Non-convergence of Newton: Theta too large
                        break; % SDNewtonLoop
                    end
                end
                NewtonIncrementOld = NewtonIncrement;

                % Update solution: Z4 <-- Z4 + DZ4;                
                Z4 = Z4 + DZ4;

                % Check error in Newton iterations
                NewtonDone = ( NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );                
                if ( NewtonDone )
                    break; % SDNewtonLoop
                end
            end

            if ( ~NewtonDone ) 
                H = Fac*H;
                Reject = true;
                SkipJac = true;
                SkipLU = false;
                continue; % Tloop;
            end
        end % End of simplified SDIRK Newton iterations

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Error estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( OPTIONS.SdirkError ) 
            DZ4 = zeros(NVAR,1);
            % L3A = 5
            if ( Coefficient.Method == 5 )
                DZ4(1:NVAR) = H*rkF(0)*F0*(1:NVAR);
                if ( rkF(1) ~= 0.0 ) 
                    DZ4 = DZ4 + rkF(1)*Z1;
                end
                if ( rkF(2) ~= 0.0 )
                    DZ4 = DZ4 + rkF(2)*Z2;
                end
                if ( rkF(3) ~= 0.0 )
                    DZ4 = DZ4 + rkF(3)*Z3;
                end
                TMP = Y + Z4;
                DZ1 = OdeFunction(T+H,TMP);
                DZ4 = DZ4 + H*rkBgam(4)*DZ1;
            else
                % DZ4(1:N) = rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 0 Z4
                if ( rkD(1) ~= 0.0 ) 
                    DZ4 = DZ4 + rkD(1)*Z1;
                end
                if ( rkD(2) ~= 0.0 )
                    DZ4 = DZ4 + rkD(2)*Z2;
                end
                if ( rkD(3) ~= 0.0 ) 
                    DZ4 = DZ4 + rkD(3)*Z3;
                end
                DZ4 = DZ4 - Z4;
            end
            
            % Perform error norm
            Err = sum((DZ4.*SCAL).^2,1);    
            Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );            
        else
            [ Err, ISTATUS ] = RK_ErrorEstimate( NVAR, H, Y, T, Z1, Z2, Z3, ...
                SCAL, FirstStep, Reject, OdeFunction, e1, ISING, Coefficient, ISTATUS );
        end
        
        % Computation of new step size Hnew
        Fac = Err^( -1/Coefficient.ELO )*min( OPTIONS.FacSafe, ( 1.0 + 2.0*OPTIONS.NewtonMaxIt ) ... 
            /( NewtonIter + 2*OPTIONS.NewtonMaxIt ) );
        Fac = min( OPTIONS.FacMax, max( OPTIONS.FacMin, Fac ) );
        Hnew = Fac*H;       

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Accept/reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ( Err < 1.0 ) % Step is accepted
            if( adjStackFlag == true )
                Zstage(1:NVAR) = Z1;
                Zstage(NVAR+1:2*NVAR) = Z2;
                Zstage(2*NVAR+1:3*NVAR) = Z3;
                stack_ptr = rkPush( NVAR, T, H, Y, Zstage, NewIt, OPTIONS, stack_ptr );                
            end
            
            if ( adjQuadFlag )
                TMP = Y + Z1;
                F = OPTIONS.QFun(T+rkC(1)*H,TMP);
                quadrature = quadrature + H*rkB(1)*F;
                
                TMP = Y + Z2;
                F = OPTIONS.QFun(T+rkC(2)*H,TMP);
                quadrature = quadrature + H*rkB(2)*F;
                
                TMP = Y + Z3;
                F = OPTIONS.QFun(T+rkC(3)*H,TMP);
                quadrature = quadrature + H*rkB(3)*F;
            end
            
            FirstStep = false;
            ISTATUS.Nacc = ISTATUS.Nacc + 1;
            if ( OPTIONS.Gustafsson ) 
                % Predictive controller of Gustafsson
                if ( ISTATUS.Nacc > 1 )
                    FacGus = OPTIONS.FacSafe*(H/Hacc)*(Err*Err/ErrOld)^(-0.25);
                    FacGus = min( OPTIONS.FacMax, max( OPTIONS.FacMin, FacGus ) );
                    Fac = min( Fac, FacGus);
                    Hnew = Fac*H;
                end
                Hacc = H;
                ErrOld = max( 1.0d-2, Err );
            end
            Hold = H;
            
            % Update time
            T = T + H;
            Tout(TYindex,1) = T;

            % Update solution: Y = Y + sum( d_i Z_i )
            if ( rkD(1) ~= 0.0 )
                Y = Y + rkD(1)*Z1;
            end
            if ( rkD(2) ~= 0.0 )
                Y = Y + rkD(2)*Z2;
            end
            if ( rkD(3) ~= 0.0 )
                Y = Y + rkD(3)*Z3;
            end        
            Yout(:,TYindex) = Y;
            
            % Update time/solution index
            TYindex = TYindex + 1;
            
            % Construct the solution quadratic interpolant Q(c_i) = Z_i. i=1:3
            if ( OPTIONS.StartNewton )
                [ Z1, Z2, Z3, CONT ] = RK_Interpolate( 'make', NVAR, H, Hold, Z1, Z2, Z3, CONT );
            end
            SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
            
            RSTATUS.Ntexit = T;
            RSTATUS.Nhnew = Hnew;
            RSTATUS.Nhacc = H;
            Hnew = Tdirection*min( max( abs(Hnew), OPTIONS.Hmin ), OPTIONS.Hmax );
            if ( Reject )
                Hnew = Tdirection*min( abs(Hnew), abs(H) );
            end
            Reject = false;
            if ( ( T + Hnew/OPTIONS.Qmin - Tfinal )*Tdirection >= 0.0 )
                H = Tfinal - T;
            else
                Hratio = Hnew/H;
                % Reuse the LU decomposition
                SkipLU = ( Theta <= OPTIONS.ThetaMin ) && ( Hratio >= OPTIONS.Qmin ) ...
                    && ( Hratio <= OPTIONS.Qmax );
                if ( ~SkipLU ) 
                    H = Hnew;
                end
            end
            % If convergence is fast enough, do not update Jacobian
            % SkipJac = ( Theta <= ThetaMin )
            SkipJac = false;

             % for debugging
             if ( OPTIONS.displaySteps == true )
                str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                disp(str);
             end
            
        else % Step is rejected
            if ( FirstStep || Reject ) 
                H = OPTIONS.FacRej*H;
            else
                H = Hnew;
            end
            Reject = true;
            SkipJac = true; % Skip if rejected - Jac is independent of H
            SkipLU = false;
            if ( ISTATUS.Nacc >= 1.0 )
                ISTATUS.Nrej = ISTATUS.Nrej + 1;
            end
            
            % for debugging
            if ( OPTIONS.displaySteps == true );
                str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                disp(str);
            end
        end
    end % Tloop
    
    % Successful exit
    Ierr = 1;
        
    % Decallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        Tout(TYindex:OPTIONS.Max_no_steps) = [];
        Yout(:,TYindex:OPTIONS.Max_no_steps) = [];       
    else
        Tout = T;
        Yout = Y;
    end
    Yout = Yout';

return;

%--------------------------------------------------------------------------

% (START) SUBROUTINE: Runge Kutta Interpolate
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Contructs or evaluates a quadratic polynomial that interpolates the Z
% solution at current step and provides the starting values for the next
% step.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ Z1, Z2, Z3, CONT ] = RK_Interpolate( action, NVAR, H, Hold, Z1, Z2, Z3, CONT )
    global rkC

    switch (action)
        case ('make')
            % Construct the solution quadratic interpolant Q(c_i) = Z_i, i=i:3
            den = ( rkC(3)-rkC(2) )*( rkC(2)-rkC(1) )*( rkC(1)-rkC(3) );
            for i=1:NVAR
                CONT(i,1) = ( -rkC(3)^2*rkC(2)*Z1(i) + Z3(i)*rkC(2)*rkC(1)^2 ...
                    + rkC(2)^2*rkC(3)*Z1(i) - rkC(2)^2*rkC(1)*Z3(i) ...
                    + rkC(3)^2*rkC(1)*Z2(i) - Z2(i)*rkC(3)*rkC(1)^2 ) ...
                    /den - Z3(i);
                CONT(i,2) = - ( rkC(1)^2*( Z3(i)-Z2(i) ) + rkC(2)^2*(Z1(i) ...
                    - Z3(i) ) + rkC(3)^2*( Z2(i)-Z1(i)) )/den;
                CONT(i,3) = ( rkC(1)*( Z3(i) - Z2(i) ) + rkC(2)*( Z1(i) - Z3(i) ) ...
                    + rkC(3)*( Z2(i) - Z1(i) ) )/den;                
            end
        case ('eval')
            % Evaluate quadratic polynomial 
            r = H/Hold;
            x1 = 1.0 + rkC(1)*r;
            x2 = 1.0 + rkC(2)*r;
            x3 = 1.0 + rkC(3)*r;
            for i=1:NVAR
                Z1(i) = CONT(i,1) + x1*( CONT(i,2) + x1*CONT(i,3) );
                Z2(i) = CONT(i,1) + x2*( CONT(i,2) + x2*CONT(i,3) );
                Z3(i) = CONT(i,1) + x3*( CONT(i,2) + x3*CONT(i,3) );
            end
    end

return; % (END) SUBROUTINE: RUNGE Kutta Interpolate

% (START) SUBROUTINE: Runge Kutta Prepare Right Hand Side
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Prepare the right-hand side for Newton iterations
% R = Z - hA*F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ R1, R2, R3, F0, L3A_flag ] = RK_PrepareRHS( T, H, Y, Z1, Z2, Z3, OdeFunction, Coefficient )
    global rkA rkC

    R1 = Z1;
    R2 = Z2;
    R3 = Z3;
    
    F0 = 'NULL';
    L3A_flag = 0;

    % L3A = 5;
    if ( Coefficient.Method == 5 )
        F0 = OdeFunction(T,Y);
        R1 = R1 - H*rkA(1,0)*F0;
        R2 = R2 - H*rkA(2,0)*F0;
        R3 = R3 - H*rkA(3,0)*F0;
        L3A_flag = 1;
    end
    
    TMP = Y + Z1;
    F = OdeFunction(T+rkC(1)*H,TMP);
    R1 = R1 - H*rkA(1,1)*F;
    R2 = R2 - H*rkA(2,1)*F;
    R3 = R3 - H*rkA(3,1)*F;
    
    TMP = Y + Z2;
    F = OdeFunction(T+rkC(2)*H,TMP);
    R1 = R1 - H*rkA(1,2)*F;
    R2 = R2 - H*rkA(2,2)*F;
    R3 = R3 - H*rkA(3,2)*F;
    
    TMP = Y + Z3;
    F = OdeFunction(T+rkC(3)*H,TMP);
    R1 = R1 - H*rkA(1,3)*F;
    R2 = R2 - H*rkA(2,3)*F;
    R3 = R3 - H*rkA(3,3)*F;
    
return; % (END) SUBROUTINE: Runge Kutaa Prepare Right Hand SIde


% (START) SUBROUTINE: Runge Kutaa Decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the matrices E1 and E2 and their decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ e1, e2, ISING, ISTATUS ] = RK_Decomp( NVAR, H, JAC, ISING, ISTATUS )
    global rkGamma rkAlpha rkBeta

    Gamma = rkGamma/H;
    Alpha = rkAlpha/H;
    Beta = rkBeta/H;
    
%~~~~~~~~~~~~~~~~~~~~TERRIBLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    % [ ISING, e1 ] = lss_decomp( NVAR, Gamma, JAC );
    e1 = -JAC + eye(NVAR,NVAR)*Gamma;  
    condNum = rcond(e1);
    if ( condNum <= eps )
        ISING = true; % ising -> 1
    else
        ISING = false; % ising -> 0
    end    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    if ( ISING ~= 0 )
        ISTATUS.Ndec = ISTATUS.Ndec;
        e2 = e1;
        return;
    end

%~~~~~~~~~~~~~~~~TERRIBLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    %[ ISING, e2 ] = lss_decomp_cmp( NVAR, Alpha, Beta, JAC );
    e2 = zeros(NVAR,NVAR);
    
    for m=1:NVAR
        for n=1:NVAR
            % e2(n,m) = -JAC(n,m) + 0*1i; % <----- confirm this matlab suggestion about imaginary representation
            e2(n,m) = -JAC(n,m); % <----- confirm this matlab suggestion about imaginary representation
        end
        cmplx = Alpha + Beta*1i; % <----- confirm this matlab suggestion about imaginary representation
        e2(m,m) = e2(m,m) + cmplx;
    end
    condNum = rcond(e2);
    if ( condNum <= eps )
        ISING = true; % ISING -> 1
    else
        ISING = false; % ISING -> 0
    end    
    ISTATUS.Ndec = ISTATUS.Ndec + 1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
return; % (END) SUBROUTINE: Runge Kutta Decomposition


% (START) SUBROUTINE: Runge Kutta Solve
function [ R1, R2, R3, ISTATUS ] = RK_Solve( NVAR, H, R1, R2, R3, e1, e2, ISING, ISTATUS )

    global rkT rkTinvAinv
    
    x(1,:) = R1./H;
    x(2,:) = R2./H;
    x(3,:) = R3./H;
    
    R = rkTinvAinv*x;
    
    R(1,:) = e1\transpose(R(1,:));

    rhs = transpose(R(2,:)) + transpose(R(3,:)).*1i;  
    
    solution = e2\rhs;
    R(2,:) = real(solution);
    R(3,:) = imag(solution);

    R = rkT*R;
    R1 = transpose(R(1,:));
    R2 = transpose(R(2,:));
    R3 = transpose(R(3,:));

    ISTATUS.Nsol = ISTATUS.Nsol + 2;

return; % (END) SUBROUTINE: Runge Kutta Solve


% (START) SUBROUTINE: Runge Kutta Error Estimate
function [ Err, ISTATUS ] = RK_ErrorEstimate( NVAR, H, Y, T, Z1, Z2, Z3, ...
    SCAL, FirstStep, Reject, OdeFunction, e1, ISING, Coefficient, ISTATUS )

    global rkE
    
    TMP = zeros(NVAR,1);
    F2 = zeros(NVAR,1);

    Transp = false;
    HEE1 = rkE(2)/H;
    HEE2 = rkE(3)/H;
    HEE3  = rkE(4)/H;
    
    F0 = OdeFunction( T, Y );
    ISTATUS.Nfun = ISTATUS.Nfun + 1;
    
    for i=1:NVAR
        F2(i) = HEE1*Z1(i) + HEE2*Z2(i) + HEE3*Z3(i);
        TMP(i) = rkE(1)*F0(i) + F2(i);
    end

    TMP = e1\TMP;
    ISTATUS.Nsol = ISTATUS.Nsol + 1;
    % R1A=2; GAU=4; L3A=5
    if ( ( Coefficient.Method == 2) || ( Coefficient.Method == 4) || ( Coefficient.Method == 5) )
        TMP = e1\TMP;
        ISTATUS.Nsol = ISTATUS.Nsol + 1;
    end
    if ( Coefficient.Method ==4 ) % Is this a wasted if statement? Deleting this does not seem to change the solution
        TMP = e1\TMP; 
        ISTATUS.Nsol = ISTATUS.Nsol + 1;
    end
    Err = errorNorm( NVAR, TMP, SCAL );
    
    if ( Err < 1.0 )
        return;
    end
    
    if ( FirstStep || Reject )
        for i=1:NVAR
            TMP(i) = TMP(i) + Y(i);
        end
        F1 = OdeFunction(T,TMP);
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        for i=1:NVAR
            TMP(i) = F1(i)+F2(i);
        end
        
        TMP = e1\TMP;
        ISTATUS.Nsol = ISTATUS.Nsol + 1;
        Err = errorNorm( NVAR, TMP, SCAL );
    end
    
return; % (END) SUBROUTINE: Runge Kutta Error Estimate

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