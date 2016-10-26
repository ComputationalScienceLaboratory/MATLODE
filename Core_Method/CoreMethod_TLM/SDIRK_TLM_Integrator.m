%% SDIRK_TLM_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr ] = SDIRK_TLM_Integrator( OdeFunction, Tspan, Y0, OPTIONS, Coefficient)
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
% Singly Diagonally Implicit Runge-Kutta tangent linear core method.
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
function [ Tout, Yout, Y_TLM, ISTATUS, RSTATUS, Ierr ] = SDIRK_TLM_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkC rkD rkE
    global rkAlpha rkGamma rkTheta

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end 

    % Get Problem Size
    NVAR = max(size(Y));
    NTLM = max(size(OPTIONS.Y_TLM));

    Z_TLM = zeros(NVAR,NTLM);
    
    Tinitial = Tspan(1);
    Tfinal = Tspan(2);
    
    Roundoff = 1.11022302462515654E-016;
    
    Yout = zeros(NVAR,OPTIONS.Max_no_steps);
    Tout = zeros(OPTIONS.Max_no_steps,1);
    TYindex = 1;
    
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    T = Tinitial;
    Tdirection = sign( Tfinal-Tinitial );
    H = max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) );
    if ( abs(H) <= 10.0*Roundoff ) 
        H = 1.0d-6;
    end
    H = min( abs(H), OPTIONS.Hmax );
    H = sign(Tdirection)*H;
    SkipLU = false;
    SkipJac = false;
    Reject = false;
    FirstStep = true;
    Transp = false;
    
    % Determine scaling factor for integration
    SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Control Flags
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GOTO_Tloop_Flag = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( (Tfinal-T)*Tdirection - Roundoff > 0.0 ) %Tloop
        % Compute E = 1/(h*gamma)-Jac and its LU decomposition
        if ( ~SkipLU )  % This time around skip the Jac update and LU
            ConsecutiveSng = 0;
            ISING = 1;
            while ( ISING ~= 0.0 ) % Hloop1
                HGammaInv = 1.0/(H*rkGamma);
                
                % Compute the Jacobian
                if ( ~SkipJac )
                    fjac = lss_jac( T,Y,OPTIONS.Jacobian );
                    ISTATUS.Njac = ISTATUS.Njac + 1;
                end
                
                [ ISING, e ] = lss_decomp( NVAR, HGammaInv, fjac );
                ISTATUS.Ndec = ISTATUS.Ndec + 1;
                
               if ( ISING ~= 0.0 )
                   str = [ 'MATRIX IS SINGULAR , ISING=', num2str(ISING), ...
                       ';     T=', num2str(T), ';     H=', num2str(H) ];
                   disp(str);
                   ISTATUS.Nsng = ISTATUS.Nsng + 1;
                   ConsecutiveSng = ConsecutiveSng + 1;
                   if ( ConsecutiveSng >= 6 )
                       return; % Failure
                   end
                   H = 0.5*H;
                   SkipJac = true;
                   SkipLU = false;
                   Reject = true;
               end
            end % Hloop1
            
            if ( ISING ~= 0.0 )
                error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;  H= %f', T, H );
            end
        end
        
        if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum bound \n T= %f;  H= %f', T, H );
        end
        if ( ( T+0.1*H == T ) || ( abs(H) <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f; H= %f', T, H );
        end
      
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Simplified Newton iterations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for istage=1:Coefficient.NStage % Stages
                        
             % Starting values for Newton iterations
             Z(:,istage) = zeros(NVAR,1);

            % Prepare the loop-independent part of the right-hand side
            G = zeros(NVAR,1);
            if ( istage > 1.0 )
                for j=1:istage-1
                    % G(:_ = sum_j Theta(i,j)*Zj(:) = H*sum_j A(i,j)*Fun(Zj(:))
                    G = G + rkTheta(istage,j)*Z(:,j);
                    % Zi(:) = sum_j Alpha(i,j)*Zj(:)
                    if ( OPTIONS.StartNewton )
                        Z(:,istage) = Z(:,istage) + rkAlpha(istage,j)*Z(:,j);
                    end
                end
            end

            % initializations for Newton iteration
            NewtonDone = false;
            Fac = 0.5;  % Step reduction factor if too many iterations

            for NewtonIter=1:OPTIONS.NewtonMaxIt % NewtonLoop
                
                % Prepare the loop-dependent part of the right-hand side
                TMP = Y + Z(:,istage);
                DZ = OdeFunction( T+rkC(istage)*H, TMP );
                ISTATUS.Nfun = ISTATUS.Nfun + 1;
                % DZ(1:N) = G(1:N) - Z(1:N. istage) + (H*rkGamma*DZ(1:N)
                DZ = DZ*H*rkGamma;
                DZ = DZ - Z(:,istage);
                DZ = DZ + G;

                % Solve the linear system
                HGammaInv = 1.0/(H*rkGamma);
                DZ = DZ*HGammaInv;
                DZ = e\DZ;
                ISTATUS.Nsol = ISTATUS.Nsol + 1;

                % Check convergence of Newton iterations
                NewtonIncrement = errorNorm( NVAR, DZ, SCAL );
                if ( NewtonIter == 1.0 )
                    Theta = abs(OPTIONS.ThetaMin);
                    NewtonRate = 2.0;
                else
                    Theta = NewtonIncrement/NewtonIncrementOld;
                    if ( Theta < 0.99 )
                        NewtonRate = Theta/(1.0-Theta);
                        % Predict error at the end of Newton process
                        NewtonPredictedErr = NewtonIncrement*Theta^( OPTIONS.NewtonMaxIt - ...
                            NewtonIter )/( 1.0 - Theta );
                        if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                            % Non-convergence of Newton: predicted error too large
                            Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS.NewtonTol );
                            Fac = 0.8*Qnewton^(-1.0/( 1 + OPTIONS.NewtonMaxIt - NewtonIter ) );
                            break; % NewtonLoop (confirm this)
                        end
                    else % Non-convergence of Newton: Theta too large
                        break; % NewtonLoop
                    end
                end
                NewtonIncrementOld = NewtonIncrement;
                % Update solution: Z(:) <-- Z(:)+DZ(:)
                Z(:,istage) = Z(:,istage) + DZ(:);

                % Check error in Newton iterations
                NewtonDone = ( NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
                if ( NewtonDone )
                    % Tune error in TLM variables by defining the minimal number of Newton iterations.
                    saveNiter = NewtonIter + 1;
                    break; % NewtonLoop (confirm this) <----------------------------------------
                end
            end % NewtonLoop

            if ( ~NewtonDone )
                H = Fac*H;
                Reject = true;
                SkipJac = true;
                SkipLU = false;
                GOTO_Tloop_Flag = true;
                break; % Tloop (confirm this: PROBLEM HERE: still in stage loop )
            end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ---> Solve for TLM variables
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Direct solution for TLM solution
            if ( OPTIONS.DirectTLM )
                TMP = Y + Z(:, istage);
                SkipJac = false;
                ConsecutiveSng = 0;
                ISING = 1;
                while ( ISING ~= 0 ) % Hloop2
                    HGammaInv = 1.0/( H*rkGamma );
                    %Compute the Jacobian
                    if ( ~SkipJac )
                        fjac1 = OPTIONS.Jacobian(T+rkC(istage)*H,TMP);
                        ISTATUS.Njac = ISTATUS.Njac + 1;
                    end
                    [ ISING, e_tlm ] = lss_decomp_tlm( NVAR, HGammaInv, ISING, fjac1 );
                    ISTATUS.Ndec = ISTATUS.Ndec + 1;
                    if ( ISING ~= 0 )
                        str = [ 'MATRIX IS SINGULAR, ISING= ', num2str(ISING), ...
                            ';     T=', num2str(T), ';     H=', num2str(H) ];
                        disp(str);
                        ISTATUS.Nsng = ISTATUS.Nsng + 1;
                        ConsecutiveSng = ConsecutiveSng + 1;
                        if ( ConsecutiveSng >= 6 )
                            return;
                        end
                        H = 0.5*H;
                        SkipJac = true;
                        SkipLU = false;
                        Reject = true;
                    end
                end

                if ( ISING ~= 0 )
                    continue; % Tloop (confirm this: PROBLEM HERE: still in stage loop )
                end

                for itlm=1:NTLM
                    G = OPTIONS.Y_TLM(:,itlm);
                    if ( istage > 1 )
                        % Gj(:) = sum_j Theta(i,j)*Zj_TLM(:)
                        %       = H*sum_j A(i,j)*Jac(Zj(:))*(Yj)TLM+Zj_TLM)
                        for j=1:istage-1
                            G = G + rkTheta(istage,j)*Z_TLM(:,j,itlm);
                        end
                    end

                    % Solve the linear system
                    HGammaInv = 1.0/(H*rkGamma);
                    G = G*HGammaInv;
                    G = e_tlm\G;
                    ISTATUS.Nsol = ISTATUS.Nsol + 1;
                    Z_TLM(:, istage, itlm) = G - OPTIONS.Y_TLM(:,itlm);
                end

            else  % DirTLM
                % Jacobian of the current stage solution
                TMP = Y + Z(:, istage);
                fjac1 = OPTIONS.Jacobian( T+rkC(istage)*H, TMP );
                ISTATUS.Njac = ISTATUS.Njac + 1;

                % Simplified Newton iterations for TLM variables
                for itlm=1:NTLM % TlmLoop
                    NewtonRate = max( NewtonRate, Roundoff )^0.8;

                    % Starting values for Newton iterations
                    Z_TLM(:,istage,itlm) = zeros(NVAR,1);

                    % Prepare the loop-independent part of the right-hand side
                    DZ = fjac1*OPTIONS.Y_TLM(:,itlm);
                    G = H*rkGamma.*DZ;
                    if ( istage > 1 )
                        % Gj(:) = sum_j Theta(i,j)*Zj_TLM(:)
                        %       = H*sum_j A(i,j)*Jac(Zj(:))*(Yj_TLM+Zj_TLM)\
                        for j=1:istage-1
                            G = G + rkTheta(istage,j)*Z_TLM(:,j,itlm);
                        end
                    end

                    % Initializations for Newton iteration
                    if ( OPTIONS.TLMNewtonEst )
                        NewtonDone = false;
                        Fac = 0.5; % Step reduction factor if too many iterations

%                         SCAL_TLM = fatOde_ErrorScale( NVAR, ITOL, OPTIONS.AbsTol_TLM(:,itlm), ...
%                             OPTIONS.RelTol_TLM(:,itlm), OPTIONS.Y_TLM(:,itlm) );
                        % Determine scaling factor for integration
                        SCAL_TLM = 1.0 ./ ( OPTIONS.AbsTol_TLM(:,itlm) + OPTIONS.RelTol_TLM(:,itlm) .* abs(OPTIONS.Y_TLM(:,itlm)) );
                    end

                    for NewtonIterTLM=1:OPTIONS.NewtonMaxIt % NewtonLoopTLM

                        % Prepare the loop-dependent part of the right-hand side
                        DZ = fjac1*Z_TLM(:,istage,itlm);
%                         DZ_temp = DZ;
%                         DZ(1:NVAR) = ( H*rkGamma )*DZ(1:NVAR) + G(1:NVAR) - ...
%                             Z_TLM(1:NVAR,istage,itlm);
                        DZ = H*rkGamma.*DZ + G - Z_TLM(:,istage,itlm);
                        

                        % Solve the linear system
                        HGammaInv = 1.0/(H*rkGamma);
                        DZ = HGammaInv.*DZ;
                        DZ = e\DZ;
                        ISTATUS.Nsol = ISTATUS.Nsol + 1;

                        if ( OPTIONS.TLMNewtonEst )
                            % Check convergence of Newton iterations
                            NewtonIncrement = errorNorm( NVAR, DZ, SCAL_TLM );
                            if ( NewtonIterTLM <= 1 )
                                ThetaTLM = abs( OPTIONS.ThetaMin );
                                NewtonRate = 2.0;
                            else
                                ThetaTLM = NewtonIncrement/NewtonIncrementOld;
                                if ( ThetaTLM < 0.99 )
                                    NewtonRate = ThetaTLM/( 1.0 - ThetaTLM );
                                    % Predict the error ar the end of Newton process
                                    NewtonPredictedErr = NewtonIncrement ...
                                        *ThetaTLM^(OPTIONS.NewtonMaxIt-NewtonIterTLM) ...
                                        /(1.0 - ThetaTLM);
                                    if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                                        % Non-convergence of Newton: predicted error tool large
                                        Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS>NewtonTol );
                                        Fac = 0.8*Qnewton^(-1.0/( 1 + OPTIONS.NewtonMaxIt ...
                                            - NewtonIterTLM ) );
                                        break; % NewtonLoopTLM (confirm this)
                                    end
                                else % Non-convergence of Newton: Theta too large
                                    break; % NewtonLoopTLM (confirm this)
                                end
                            end
                            NewtonIncrementOld = NewtonIncrement;
                        end % TLMNewtonEst

                        % Update solution: Z_TLM(:) <-- Z_TLM(:) + DZ(:)
                        Z_TLM(:,istage,itlm) = Z_TLM(:,istage,itlm) + DZ(:);

                        % Check error in Newton iterations
                        if ( OPTIONS.TLMNewtonEst )
                            NewtonDone = (NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
                            if ( NewtonDone )
                                break; % NewtonLoopTLM (confirm this)
                            end
                        else
                            % Minimum number of iterations same as FWD iterations 
                            if ( NewtonIterTLM >= saveNiter )
                                break; % NewtonLoopTLM (confirm this)
                            end
                        end
                    end % NewtonLoopTLM

                    if ( OPTIONS.TLMNewtonEst && (~NewtonDone) )
                        H = Fac*H;
                        Reject  = true;
                        SkipJac = true;
                        SkipLU  = false;
                        GOTO_Tloop_Flag = true;
                        break; % Tloop (**confirm this)
                    end
                end % TlmLoop
            end % DirTLM
        end % Stages
        if ( GOTO_Tloop_Flag == true )
            GOTO_Tloop_Flag = false;
            continue; %Tloop
        end
                        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Error Estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        Yerr = zeros(NVAR,1);
        for i=1:Coefficient.NStage
            if ( rkE(i) ~= 0.0 )
                Yerr = Yerr + rkE(i)*Z(:,i);
            end
        end

        % Solve the linear system
        HGammaInv = 1.0/( H*rkGamma );
        Yerr = Yerr*HGammaInv;
        Yerr = e\Yerr;
        ISTATUS.Nsol = ISTATUS.Nsol + 1;

        Err = errorNorm( NVAR, Yerr, SCAL );

        if (OPTIONS.TLMTruncErr)
            Yerr_TLM = zeros(NTLM,1);
            for itlm=1:NTLM
                for j=1:Coefficient.NStage
                    if ( rkE(j) ~= 0 )
                        Yerr_TLM(:,itlm) = Yerr_TLM(:,itlm) + rkE(j)*Z_TLM(:,j,itlm);
                    end
                    Yerr_TLM(:,itlm) = HGammaInv.*Yerr_TLM(:,itlm);
                    Yerr_TLM = e\Yerr_TLM;
                    ISTATUS.Nsol = ISTATUS.Nsol + 1;
                end
                Err = SDIRK_ErrorNorm_TLM(NVAR, NTLM, Yerr_TLM, OPTIONS.AbsTol_TLM, OPTIONS.RelTol_TLM, Err );
            end
        end
        
        % Computation of new step size Hnew
        Fac = OPTIONS.FacSafe*Err^(-1.0/Coefficient.ELO);
        Fac = max( OPTIONS.FacMin, min(OPTIONS.FacMax, Fac ) );
        Hnew = H*Fac;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Accept/Reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( Err < 1.0 ) % accept
            FirstStep = false;
            ISTATUS.Nacc = ISTATUS.Nacc + 1;

            % Update time
            T = T + H;
            
            % Update solution
            for i=1:Coefficient.NStage
                if ( rkD(i) ~= 0.0 )
                    Y = Y + rkD(i)*Z(:,i);
                    for itlm=1:NTLM
                        OPTIONS.Y_TLM(:,itlm) = OPTIONS.Y_TLM(:,itlm) + rkD(i)*Z_TLM(:,i,itlm);
                    end
                end
            end
            
            % Store checkpoint values
            if ( OPTIONS.storeCheckpoint == true )
                Tout(TYindex,1) = T;            
                Yout(:,TYindex) = Y;
                TYindex = TYindex + 1;
            end
            
            % Determine scaling factor for integration
            SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );

            % Next time step
            Hnew = Tdirection*min( abs(Hnew), OPTIONS.Hmax );
            
            % Last T and H
            RSTATUS.Ntexit = T;
            RSTATUS.Nhexit = H;
            RSTATUS.Nhexit = Hnew;
            
            % No step increase after rejection
            if ( Reject ) 
                Hnew = Tdirection*min( abs(Hnew), abs(H) );
            end
            Reject = false;
            
            if ( ( T + Hnew/OPTIONS.Qmin - Tfinal )*Tdirection > 0.0 )
                H = Tfinal-T;
            else
                Hratio = Hnew/H;
                % if step not changed too much keep Jacobian and reuse LU
                % SkipLU = ( ( Theta <= ThetaMin ) && ( Hratio >= Qmin ) && ( Hratio <= Qmax ) );
                % For TLM: do not skip LU (decrease TLM error)
                SkipLU = false;
                if ( ~SkipLU )
                    H = Hnew;
                end
            end
            
            % If convergence is fast enough, do not update Jacobian
            SkipJac = false;
            
            % for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Accpeted step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);            
            end
            
        else % reject
            if ( FirstStep || Reject )
                H = OPTIONS.FacRej*H;
            else
                H = Hnew;
            end
            Reject  = true;
            SkipJac = true;
            SkipLU  = false;
            if ( ISTATUS.Nacc >= 1 )
                ISTATUS.Nrej = ISTATUS.Nrej + 1;
            end
            
            %for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Rejected step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        end % accept
    end % Tloop
    
    Y_TLM = OPTIONS.Y_TLM;
    
    % Successful return
    Ierr = 1;
    
    if ( OPTIONS.storeCheckpoint == true )
        Tout(TYindex:OPTIONS.Max_no_steps,:) = [];
        Yout(:,TYindex:OPTIONS.Max_no_steps) = [];
    else
        Tout = T;
        Yout = Y;
    end
       
    % Temporary Fix
    Yout = transpose(Yout);
    
return;

%--------------------------------------------------------------------------

% (START) SUBROUNTINE: Singly Diagonally Impliit Runge Kutta Error Norm Tangent Linear Model
function [ FWD_Err ] = SDIRK_ErrorNorm_TLM( NVAR, NTLM, Y_TLM, AbsTol_TLM, RelTol_TLM, FWD_Err )
    for itlm=1:NTLM
        SCAL_TLM = fatOde_ErrorScale( NVAR, 1, AbsTol_TLM(:,itlm), RelTol_TLM(:,itlm), ...
            Y_TLM(:,itlm) );
        Err = errorNorm( NVAR, Y_TLM(:,itlm), SCAL_TLM );
        FWD_Err = max( FWD_Err, Err );
    end

return; % (END) SUBROUTINE: Singly Diagonally Implicit Runge Kutta Error Norm Tangent Linear Model


