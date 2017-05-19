%% RK_TLM_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr ] = RK_TLM_Integrator( OdeFunction, Tspan, Y0, OPTIONS, Coefficient)
%
%% Input Parameters
% |OdeFunction|: ODE function function handle
%
% |Tspan|: Time interval
%
% |Y0|: Initial state vector
%
% |OPTIONS|: Option struct
%
% |Coefficients|: Constant coefficients associated with method
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
%% Description
% Implicit Runge-Kutta tangent linear core method.
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
function [ Tout, Yout, Y_TLM, ISTATUS, RSTATUS, Ierr ] = RK_TLM_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient )
    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end  

    % Get Problem Size
    NVAR = max(size(Y));
    NTLM = max(size(OPTIONS.Y_TLM));

    Tinitial = Tspan(1);
    Tfinal   = Tspan(2);
    
    Roundoff = 1.11022302462515654E-016;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkC rkD
    global rkBgam 
    global rkGamma rkTheta
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONT = zeros(NVAR,3);
    
%    AbsTol_TLM = AbsTol; % <--- Temporary Fix
%    RelTol_TLM = RelTol; % <--- Temporary Fix
    
    Yout = zeros(NVAR, OPTIONS.Max_no_steps);
    Tout = zeros(OPTIONS.Max_no_steps,1);
    TYindex = 1;
    
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;

    saveNiter = 0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Required
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    if ( isempty(OPTIONS.Y_TLM) )
        error('User defined option parameter Y_TLM is required.');
    end    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%OPTIONS.SdirkError = 0;  
%OPTIONS.TLMTruncErr = 0;
%OPTIONS.Gustafsson = 0;
%OPTIONS.TLMNewtonEst = 0;

    T = Tinitial;
    Tdirection = sign( Tfinal-Tinitial );
    H = min( max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) ), OPTIONS.Hmax );
    if ( abs(H) <= 10.0*Roundoff )
        H = 1.0d-6;
    end
    H           = H*sign(Tdirection);
    Hold        = H;
    Reject      = false;
    FirstStep   = true;
    SkipJac     = false;
    SkipLU      = false;
    Transp      = false;
    ISING       = false;
    if ( (Tinitial+H*1.0001-Tfinal)*Tdirection >= 0.0 )
        H = Tfinal-T;
    end
    Nconsecutive = 0;
    
    Ierr = 0;
    
    % Determine scaling factor for integration
    SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
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
                    error('Non-convergence of Newton iterations');
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
            error('Number of steps exceeds maximum bound');
        end

        if ( 0.1*abs(H) <= abs(T)*Roundoff )
            error('Number of steps exceeds maximum bound');
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
            [ DZ1, DZ2, DZ3, F0, L3A_flag, ISTATUS ] = RK_PrepareRHS( T, H, Y, Z1, Z2, Z3, OdeFunction, Coefficient, ISTATUS );
            
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
                saveNiter = NewtonIter + 1;
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
            DZ4 = OdeFunction(T,Y);
            ISTATUS.Nfun = ISTATUS.Nfun + 1;
            
            % G = H*rkBgam(0)*F0 + rkTheta(1)*Z1 + rkTheta(2)*Z2 + rkTheta(3)*Z3
            G = zeros(NVAR,1);
            G = G + rkBgam(1)*H*DZ4;
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
%   Error estimation, forward solution
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( OPTIONS.SdirkError )
            % DZ4(1:N) =  rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 - Z4
            DZ4 = zeros(NVAR,1);
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
            Err = errorNorm( NVAR, DZ4, SCAL );
        else
            [ Err, ISTATUS ] = RK_ErrorEstimate( NVAR, H, Y, T, Z1, Z2, Z3, ...
                SCAL, FirstStep, Reject, OdeFunction, e1, ISING, Coefficient, ISTATUS );
        end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  If error small enough, compute TLM solution
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( Err < 1.0 )
            
            % Jacobian values
            TMP = Y + Z1;
            fjac1 = OPTIONS.Jacobian(T+rkC(1)*H,TMP);
            TMP = Y + Z2;
            fjac2 = OPTIONS.Jacobian(T+rkC(2)*H,TMP);
            TMP = Y + Z3;
            fjac3 = OPTIONS.Jacobian(T+rkC(3)*H,TMP);
            ISTATUS.Njac = ISTATUS.Njac + 3;           
            
            if ( OPTIONS.DirectTLM ) % TLMDIR
                % Construct the full Jacobian
                ax_big = lss_decomp_big( NVAR, H, fjac1, fjac2, fjac3 );
                ISTATUS.Ndec = ISTATUS.Ndec + 1;
 
                Z1_TLM = zeros(NVAR,NTLM);
                Z2_TLM = zeros(NVAR,NTLM);
                Z3_TLM = zeros(NVAR,NTLM);
                for itlm = 1:NTLM
                    Zbig = RK_PrepareRHS_TLMdirect( NVAR, H, OPTIONS.Y_TLM(:,itlm), fjac1, fjac2, fjac3);
                    Zbig = Zbig'; %%%%%% Temporary fix
                    Zbig = ax_big\Zbig;
                    Z1_TLM(:,itlm) = Zbig(1:NVAR);
                    Z2_TLM(:,itlm) = Zbig(NVAR+1:2*NVAR);
                    Z3_TLM(:,itlm) = Zbig(2*NVAR+1:3*NVAR);
                end                   
                ISTATUS.Nsol = ISTATUS.Nsol + NTLM;
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Loop for Newton iterations, TLM variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else % TLMDIR
                for itlm=1:NTLM %Tlm
                    % Starting values for Newton iterations
                    Z1_TLM(:,itlm) = zeros(NVAR,1);
                    Z2_TLM(:,itlm) = zeros(NVAR,1);
                    Z3_TLM(:,itlm) = zeros(NVAR,1);
                    
                    % Initializations for Newton iteration
                    if ( OPTIONS.TLMNewtonEst )
                        NewtonDone = false;
                        Fac = 0.5; % Step reduction if too many iterations
                        SCAL_TLM = fatOde_ErrorScale( NVAR, OPTIONS.ITOL, OPTIONS.AbsTol_TLM(:,itlm), ...
                            OPTIONS.RelTol_TLM(:,itlm), OPTIONS.Y_TLM(:,itlm) );
                        
                        % Determine scaling factor for integration
%                         SCAL_TLM = 1.0 ./ ( OPTIONS.AbsTol_TLM(:,itlm) + OPTIONS.RelTol_TLM(:,itlm) .* abs(OPTIONS.Y_TLM(:,itlm)) );                        
                    end
                    
                    for NewtonIterTLM=1:OPTIONS.NewtonMaxIt % NewtonLoopTLM
                        % Prepare the right-hand side
                        [ DZ1, DZ2, DZ3 ] = RK_PrepareRHS_TLM( NVAR, H, OPTIONS.Y_TLM(:,itlm), ...
                            Z1_TLM(:,itlm), Z2_TLM(:,itlm), Z3_TLM(:,itlm), ...
                            DZ1, DZ2, DZ3, fjac1, fjac2, fjac3 );
                        
                        % Solve the linear systems
                        [ DZ1, DZ2, DZ3, ISTATUS ] = RK_Solve( NVAR, H, DZ1, DZ2, DZ3, e1, e2, ISING, ISTATUS);
                        if ( OPTIONS.TLMNewtonEst )
                            % Check convergence of Newton iterations
                            NewtonIncrement = sqrt( ( RK_ErrorNorm( NVAR, SCAL_TLM, DZ1)^2 + ...
                                RK_ErrorNorm( NVAR, SCAL_TLM, DZ2 )^2 + ...
                                RK_ErrorNorm( NVAR, SCAL_TLM, DZ3 )^2 )/3.0 );
                                                
                            if ( NewtonIterTLM == 1.0 )
                                ThetaTLM = abs( OPTIONS.ThetaMin );
                                NewtonRate = 2.0;
                            else
                                ThetaTLM = NewtonIncrement/NewtonIncrementOld;
                                if ( ThetaTLM < 0.99 )
                                    NewtonRate = ThetaTLM/(1.0-ThetaTLM);
                                else % Non-convergence of Newton: ThetaTLM too large
                                    break; % NewtonLoopTLM (confirm this)
                                end
                                if ( NewtonIterTLM < OPTIONS.NewtonMaxIt )
                                    % Predict error at the end of Newton process
                                    NewtonPredictedErr = NewtonIncrement + ...
                                        ThetaTLM^(OPTIONS.NewtonMaxIt-NewtonIterTLM)/ ...
                                        (1.0-ThetaTLM);
                                    if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                                        % Non-convergence of Newton: predicted error too large
                                        Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS.NewtonTol );
                                        Fac = 0.8*Qnewton^( -1.0/(10+OPTIONS.NewtonMaxIt-NewtonIterTLM) );
                                        break; % NewtonLoopTlM
                                    end
                                end
                            end
                            NewtonIncrementOld = max( NewtonIncrement, Roundoff );
                        end % TLMNewtonEst
                        
                        % Update solution
                        %Z1 = Z1 - DZ1;
                        %Z2 = Z2 - DZ2;
                        %Z3 = Z3 - DZ3;
                        Z1_TLM(:,itlm) = Z1_TLM(:,itlm) - DZ1;
                        Z2_TLM(:,itlm) = Z2_TLM(:,itlm) - DZ2;
                        Z3_TLM(:,itlm) = Z3_TLM(:,itlm) - DZ3;
                        
                        % Check error in Newton iterations
                        if ( OPTIONS.TLMNewtonEst )
                            NewtonDone = ( NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
                            if ( NewtonDone )
                                break; % newtonLoopTLM 
                            end
                        else
                            % Minimum number of iteations same as FWD iterations
                            if ( NewtonIterTLM >= saveNiter )
                                break; % NewtonLoopTLM 
                            end
                        end
                        
                    end
                    
                    if ( OPTIONS.TLMNewtonEst && (~NewtonDone) )
                        H = Fac*H;
                        Reject = true;
                        SkipJac = true;
                        SkipLU = false;
                        continue; % check this
                    end
                end
            end
                        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error estimation, TLM solution
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if ( OPTIONS.TLMTruncErr )
                [ Err, ISTATUS ] = RK_ErrorEstimate_TLM( NVAR, NTLM, T, H, Y, OPTIONS, ...
                    Z1_TLM, Z2_TLM, Z3_TLM, Err, FirstStep, Reject, ISTATUS );
            end
        end % ( Err<1.0 )
        
        % Computation of new step size Hnew
        Fac = Err^( -1/Coefficient.ELO )*min( OPTIONS.FacSafe, ( 1.0 + 2.0*OPTIONS.NewtonMaxIt ) ... 
            /( NewtonIter + 2*OPTIONS.NewtonMaxIt ) );
        Fac = min( OPTIONS.FacMax, max( OPTIONS.FacMin, Fac ) );
        Hnew = Fac*H;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Accept/reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( Err < 1.0 ) % STEP IS ACCEPTED
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
            
            % Update TLM solution: Y <- Y + sum(d_i*Z_i_TLM)
            for itlm = 1:NTLM
                if ( rkD(1) ~= 0.0 )
                    OPTIONS.Y_TLM(:,itlm) = OPTIONS.Y_TLM(:,itlm) + rkD(1)*Z1_TLM(:,itlm);
                end
                if ( rkD(2) ~= 0.0 )
                    OPTIONS.Y_TLM(:,itlm) = OPTIONS.Y_TLM(:,itlm) + rkD(2)*Z2_TLM(:,itlm);
                end
                if ( rkD(3) ~= 0.0 )
                    OPTIONS.Y_TLM(:,itlm) = OPTIONS.Y_TLM(:,itlm) + rkD(3)*Z3_TLM(:,itlm);
                end
            end
            
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
            
            % for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Accpeted step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        else % Skip is rejected
            if ( FirstStep || Reject )
                H = OPTIONS.FacRej*H;
            else
                H = Hnew;
            end
            
            Reject = true;
            SkipJac = true;
            SkipLU = false;
            
            % Update Number of rejections
            if ( ISTATUS.Nacc >= 1.0 )
                ISTATUS.Nrej = ISTATUS.Nrej + 1;
            end
            
            %for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Rejected step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        end
    end % Tloop
    
    Tout(TYindex:OPTIONS.Max_no_steps,:) = [];
    Yout(:,TYindex:OPTIONS.Max_no_steps) = [];
    
    Yout = transpose(Yout);
    
    Y_TLM = OPTIONS.Y_TLM;
    
return;

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

return; 

function [ Err ] = RK_ErrorNorm(NVAR,SCAL,TMP)
    % Perform error norm
    Err = sum((TMP.*SCAL).^2,1);
    Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );   
return;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Prepare the right-hand side for Newton iterations
% R = Z - hA*F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ R1, R2, R3, F0, L3A_flag, ISTATUS ] = RK_PrepareRHS( T, H, Y, Z1, Z2, Z3, OdeFunction, Coefficient, ISTATUS )
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
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
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
    
    ISTATUS.Nfun = ISTATUS.Nfun + 3;
    
return;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the matrices E1 and E2 and their decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ e1, e2, ISING, ISTATUS ] = RK_Decomp( NVAR, H, JAC, ISING, ISTATUS )
    global rkGamma rkAlpha rkBeta

    Gamma = rkGamma/H;
    Alpha = rkAlpha/H;
    Beta = rkBeta/H;
    [ ISING, e1 ] = lss_decomp( NVAR, Gamma, JAC );
    
    if ( ISING ~= 0 )
        ISTATUS.Ndec = ISTATUS.Ndec;
        e2 = e1;
        return;
    end

    [ ISING, e2 ] = lss_decomp_cmp( NVAR, Alpha, Beta, JAC );
    ISTATUS.Ndec = ISTATUS.Ndec + 1;
    
return; 


function [ R1, R2, R3, ISTATUS ] = RK_Solve( NVAR, H, R1, R2, R3, e1, e2, ISING, ISTATUS )

    global rkT rkTinvAinv

    x1 = R1./H;
    x2 = R2./H;
    x3 = R3./H;
    R1 = rkTinvAinv(1,1).*x1 + rkTinvAinv(1,2).*x2 + rkTinvAinv(1,3).*x3;
    R2 = rkTinvAinv(2,1).*x1 + rkTinvAinv(2,2).*x2 + rkTinvAinv(2,3).*x3;
    R3 = rkTinvAinv(3,1).*x1 + rkTinvAinv(3,2).*x2 + rkTinvAinv(3,3).*x3;
   
    R1 = e1\R1;

    rhs = R2 + R3.*1i;
    
    solution = e2\rhs;
    R2 = real(solution);
    R3 = imag(solution);
        
    x1 = R1;
    x2 = R2;
    x3 = R3;
    R1 = rkT(1,1).*x1 + rkT(1,2).*x2 + rkT(1,3).*x3;
    R2 = rkT(2,1).*x1 + rkT(2,2).*x2 + rkT(2,3).*x3;
    R3 = rkT(3,1).*x1 + rkT(3,2).*x2 + rkT(3,3).*x3;    
    
    ISTATUS.Nsol = ISTATUS.Nsol + 2;

return; 


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    
return;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ FWD_Err, ISTATUS ] = RK_ErrorEstimate_TLM( NVAR, NTLM, T, H, Y, OPTIONS, ...
    Z1_TLM, Z2_TLM, Z3_TLM, FWD_Err, FirstStep, Reject, ISTATUS )

    global rkE
    
    Transp = false;
    HEE1 = rkE(2)/H;
    HEE2 = rkE(3)/H;
    HEE3 = rkE(4)/H;
    
    for itlm=1:NTLM
        SCAL_TLM = 1.0 ./ ( OPTIONS.AbsTol_TLM(:,itlm) + OPTIONS.RelTol_TLM(:,itlm) .* abs(OPTIONS.Y_TLM(:,itlm)) );
        
        fjac = OPTIONS.Jacobian(T,Y);
        ISTATUS.Njac = ISTATUS.Njac + 1;
        OPTIONS.Y_TLM(:,itlm) = fjac*JY_TLM;
        
        for i=1:NVAR
            TMP2(i) = HEE1*Z1_TLM(:,itlm) + HEE2*Z2_TLM(i,itlm) + ...
                HEE3*Z3_TLM(i,itlm);
            TMP(i) = rkE(1)*JY_TLM(i) + TMP2(i);
        end
        
        lss_solve( Transp, TMP );
        ISTATUS.Nsol = ISTATUS.Nsol + 1;
        
        Err = RK_ErrorNorm( NVAR, SCAL_TLM, TMP );
        
        FWD_Err = max( FWD_Err, Err );
    end
        
    
return; % (END) SUBROUTINE: Runge-Kutta Error Estimate Tangent Linear Model

% (START) SUBROUTINE: Runge Kutta Prepare Right Hand Side Tangent Linear Model
function [ R1 R2 R3 ] = RK_PrepareRHS_TLM( NVAR, H, Y_TLM, ...
    Z1_TLM, Z2_TLM, Z3_TLM, R1, R2, R3, fjac1, fjac2, fjac3 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Prepare the right-hand side for Newton iterations
%   R = Z_TLM - hA * Jac*Z_TLM
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global rkA

    R1 = Z1_TLM;
    R2 = Z2_TLM;
    R3 = Z3_TLM;
    
    TMP = Y_TLM + Z1_TLM;
    F = fjac1*TMP; 
    R1 = R1 - H*rkA(1,1)*F;
    R2 = R2 - H*rkA(2,1)*F;
    R3 = R3 - H*rkA(3,1)*F;
    
    TMP = Y_TLM + Z2_TLM;
    F = fjac2*TMP; 
    R1 = R1 - H*rkA(1,2)*F;
    R2 = R2 - H*rkA(2,2)*F;
    R3 = R3 - H*rkA(3,2)*F;
    
    TMP = Y_TLM + Z3_TLM;
    F = fjac3*TMP; 
    R1 = R1 - H*rkA(1,3)*F;
    R2 = R2 - H*rkA(2,3)*F;
    R3 = R3 - H*rkA(3,3)*F;

return; % (END) SUBROUTINE: Runge Kutta Prepare Right Hand Side Tangent Linear Model

% (START) SUBROUTINE: Runge Kutta Prepare Right Hand Side Tangent Linear Model Direct
function [ Zbig ] = RK_PrepareRHS_TLMdirect( NVAR, H, Y_TLM, fjac1, fjac2, fjac3 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Prepare the right-hand side for direct solution
%   Z = hA * Jac*Y_TLM
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global rkA

    F1 = fjac1*Y_TLM;
    F2 = fjac2*Y_TLM;
    F3 = fjac3*Y_TLM;
    % ---> Commented code for sparse big linear algebra:
    % Zbig(1,1:NVAR) = H*rkA(1,1)*F(1:NVAR);
    % Zbig(2,1:NVAR) = H*rkA(2,1)*F(1:NVAR);
    % Zbig(3,1:NVAR) = H*rkA(3,1)*F(1:NVAR);
    Zbig(1:NVAR)          = H*rkA(1,1)*F1(1:NVAR) + H*rkA(1,2)*F2(1:NVAR) + H*rkA(1,3)*F3(1:NVAR);
    Zbig(NVAR+1:2*NVAR)   = H*rkA(2,1)*F1(1:NVAR) + H*rkA(2,2)*F2(1:NVAR) + H*rkA(2,3)*F3(1:NVAR);
    Zbig(2*NVAR+1:3*NVAR) = H*rkA(3,1)*F1(1:NVAR) + H*rkA(3,2)*F2(1:NVAR) + H*rkA(3,3)*F3(1:NVAR);
    
    %F = fjac2*Y_TLM;
    % ---> Commented cade for sparse big linear algebra:
    % Zbig(1,1:NVAR) = Zbig(1,1:NVAR) + H*rkA(1,2)*F(1:NVAR);
    % Zbig(2,1:NVAR) = Zbig(2,1:NVAR) + H*rkA(2,2)*F(1:NVAR);
    % Zbig(3,1:NVAR) = Zbig(3,1:NVAR) + H*rkA(3,2)*F(1:NAR);
    %Zbig(1:NVAR)          = Zbig(1:NVAR)          + H*rkA(1,2)*F2(1:NVAR);
    %Zbig(N+1:2*NVAR)      = Zbig(N+1:2*NVAR)      + H*rkA(2,2)*F2(1:NVAR);
    %Zbig(2*NVAR+1:3*NVAR) = Zbig(2*NVAR+1:3*NVAR) + H*rkA(3,2)*F2(1:NVAR);
    
    %F = fjac3*Y_TLM;
    % ---> Commented code for sparse big linear algebra:
    % Zbig(1,1:NVAR) = Zbig(1,1:NVAR) + H*rkA(1,3)*F(1:NVAR);
    % Zbig(2,1:NVAR) = Zbig(2,1:NVAR) + H*rkA(2,3)*F(1:NVAR);
    % Zbig(3,1:NVAR) = Zbig(3,1:NVAR) + H*rkA(3,3)*F(1:NVAR);
    %Zbig(1:NVAR)          = Zbig(1:NVAR)          + H*rkA(1,3)*F3(1:NVAR);
    %Zbig(NVAR+1:2*NVAR)   = Zbig(N+1:2*NVAR)      + H*rkA(2,3)*F3(1:NVAR);
    %Zbig(2*NVAR+1:3*NVAR) = Zbig(2*NVAR+1:3*NVAR) + H*rkA(3,3)*F3(1:NVAR);

return; % (END) SUBROUNTINE: Runge Kutta Prepare Right Hand Side Tangent Linear Model Direct


