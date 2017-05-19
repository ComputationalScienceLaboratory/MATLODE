%% SDIRK_FWD_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = SDIRK_FWD_Integrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )
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
% Singly Diagonally Implicit Runge-Kutta forward core method.
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
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = SDIRK_FWD_Integrator( OdeFunction,...
    Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )

    % Get Problem Size
    NVAR = max(size(Y));

    Tinitial = Tspan(1);
    Tfinal = Tspan(2);

    Roundoff = 1.11022302462515654E-016;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkGamma
    global rkB
    global rkC
    global rkTheta
    global rkAlpha
    global rkE
    global rkD

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gmresFlag = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end

    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');

    quadrature = OPTIONS.Quadrature();
    
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
    
    SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
    
   if ( OPTIONS.storeCheckpoint == true )
        Yout = zeros(NVAR,OPTIONS.Max_no_steps);
        Tout = zeros(OPTIONS.Max_no_steps,1);
        TYindex = 1;
        Yout(:,TYindex) = Y;
        Tout(TYindex,1) = Tinitial; 
    end
    
    % Preallocate for speed
    Z = zeros(NVAR,Coefficient.NStage);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Control Flags
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GOTO_Tloop_Flag = false;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Tfinal-T )*Tdirection - Roundoff > 0.0 )
        % Compute E = 1/(h*gamma)-Jac and its LU decomposition
        if ( ~SkipLU )
            ConsecutiveSng = 0;
            ISING = 1;
            while ( ISING ~= 0.0 ) % Hloop1
                HGammaInv = 1.0/(H*rkGamma);
                
                % Compute the Jacobian
                if ( ~SkipJac )
                    if ( ~OPTIONS.MatrixFree )
                        fjac = OPTIONS.Jacobian(T,Y);
                        ISTATUS.Njac = ISTATUS.Njac + 1;
                    else
                        if( ~isempty(OPTIONS.Jacobian) )
                            fjac = @(vee)OPTIONS.Jacobian(T,Y,vee);
                        else
                            Fcn0 = OdeFunction(T,Y);
                            normy = norm(Y);
                            fjac = @(v)Mat_Free_Jac(T,Y,v,OdeFunction,Fcn0,normy);
                        end
                    end
                end
                
                % Decomposition
                if ( ~OPTIONS.MatrixFree )
                    e = -fjac + speye(NVAR,NVAR)*HGammaInv; % faster since sparse
                    condNum = rcond(full(e));
                else
                    e = @(v)( HGammaInv*v -fjac(v));
                    condNum = 1; % Place Holder
                end
                
                if ( condNum < 2.2204e-14 || gmresFlag ~= 0 ) % 2.2204e-14 = 100*eps
                    ISING = true; % ising -> 1
                else
                    ISING = false; % ising -> 0
                end                
                ISTATUS.Ndec = ISTATUS.Ndec + 1;
                
               if ( ISING ~= 0.0 )
                   str = [ 'MATRIX IS SINGULAR , ISING=', num2str(ISING), ...
                       ';     T=', num2str(T), ';     H=', num2str(H) ];
                   disp(str);
                   ISTATUS.Nsng = ISTATUS.Nsng + 1;
                   ConsecutiveSng = ConsecutiveSng + 1;
                   if ( ConsecutiveSng >= 6 )
                       Ierr = 1;
                       return; % Failure
                   end
                   H = 0.5*H;
                   SkipJac = true;
                   SkipLU = false;
                   Reject = true;
                   
                   if ( gmresFlag ~= 0 )
                       ISING = 0;
                       break;
                   end

               end
            end % Hloop1
           
            if ( ISING ~= 0 )
                error('Matrix is repeatedly singular');
            end
%         elseif(OPTIONS.MatrixFree)
%             HGammaInv = 1.0/(H*rkGamma);
%             fjac = @(v)OPTIONS.Jacobian(T,Y,v);
%             e = @(v)(v - HGammaInv*fjac(v));
         end
    
        if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum bound: %d',ISTATUS.Nstp);
        end

        if ( (T+0.1*H == T) || (abs(H) <= Roundoff) )
            error('Step size too small: T + 10*H = T or H < Roundoff');
        end
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for istage=1:Coefficient.NStage % stages
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Simplified Newton iterations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                if( ~OPTIONS.MatrixFree )
                    if ( OPTIONS.LU )
                        %%%%% LU Decomp NEW %%%%%
                        [L,U,P] = lu(e);
                        DZ = U\(L\(P*DZ));
                        %%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                        DZ = e\DZ;
                    end                                       
                else
                    [ DZ, gmresFlag, ~, iter ] = gmres(e, DZ, [], ...
                        OPTIONS.GMRES_TOL);
                    ISTATUS.Nfun = ISTATUS.Nfun + iter(1);
                    switch(gmresFlag)
                        case 1
                            warning('GMRES: iterated MAXIT times but did not converge');
                            break;
                        case 2
                            warning('GMRES: preconditioner M was ill-conditioned');
                            break;
                        case 3
                            warning('GMRES: stagnated (two consecutive iterates were the same)');
                            break;
                    end
                end
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

                % check error in Newton iterations
                NewtonDone = (NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
                if ( NewtonDone )
                    break; % NewtonLoop
                end
            end % NewtonLoop
            
            if ( gmresFlag ~= 0 )
                break;
            end

%             if ( ISTATUS.Nstp == 5 )
%                 keyboard
%             end
            
            if ( ~NewtonDone )
                H = Fac*H;
                Reject = true;
                SkipJac = true;
                SkipLU = false;
                GOTO_Tloop_Flag = true;
                break; % Tloop
            end
        end % Stages 
        
        if ( gmresFlag ~= 0 )
            continue;
        end
        
        if ( GOTO_Tloop_Flag == true )
            GOTO_Tloop_Flag = false;
            continue; % Tloop
        end;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error Estimation
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
        if(~OPTIONS.MatrixFree)
            if ( OPTIONS.LU )
                %%%%% LU Decomp NEW %%%%%
                [L,U,P] = lu(e);
                Yerr = U\(L\(P*Yerr));
                %%%%%%%%%%%%%%%%%%%%%%%%%                
            else
                Yerr = e\Yerr;
            end                        
        else
            [ Yerr, gmresFlag, ~, iter ] = gmres(e, Yerr, ...
                [], OPTIONS.GMRES_TOL);
            switch(gmresFlag)
                case 1
                    warning('GMRES: iterated MAXIT times but did not converge');
                    break;
                case 2
                    warning('GMRES: preconditioner M was ill-conditioned');
                    break;
                case 3
                    warning('GMRES: stagnated (two consecutive iterates were the same)');
                    break;
            end            
            ISTATUS.Nfun = ISTATUS.Nfun + iter(1);
        end
        ISTATUS.Nsol = ISTATUS.Nsol + 1;       

        % Calculate error norm
        Err = sum((Yerr.*SCAL).^2,1);    
        Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );

        % Computation of new step size Hnew
        Fac = OPTIONS.FacSafe*(Err)^(-1.0/Coefficient.ELO);
        Fac = max(OPTIONS.FacMin, min( OPTIONS.FacMax, Fac ));
        Hnew = H*Fac;
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Accept/Reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( Err < 1.0 )    % accepted
            FirstStep = false;
            ISTATUS.Nacc = ISTATUS.Nacc + 1;

            % Checkpoints for Adjoint Calculation
            if ( adjStackFlag == 1 )
                stack_ptr = sdirkPush(NVAR,Coefficient.NStage,T,H,Y,Z,OPTIONS.Max_no_steps,stack_ptr);
            end            
            
            % Update the results for the quadrature term
            if ( adjQuadFlag )
                for i=1:Coefficient.NStage
                    TMP = Y + Z(:,i);
                    R = OPTIONS.QFun(T+rkC(i)*H,TMP);
                    quadrature = quadrature + H*rkB(i)*R;
                end
            end
            
            % Update time
            T = T + H;
           
            % Update solution
            for i=1:Coefficient.NStage
                if ( rkD(i) ~= 0.0 )
                    Y = Y + rkD(i)*Z(:,i);
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

            % Next Time step
            Hnew = Tdirection*min(abs(Hnew),OPTIONS.Hmax);
            
            % Last T and H
            RSTATUS.Ntexit = T;
            RSTATUS.Nhexit = H;
            RSTATUS.Nhnew = Hnew;
            
            % No step increase after a rejection
            if ( Reject )
                Hnew = Tdirection*min(abs(Hnew),abs(H));
            end
            Reject = false;
            if ( (T+Hnew/OPTIONS.Qmin-Tfinal)*Tdirection > 0.0 )
                H = Tfinal-T;
            else
                Hratio = Hnew/H;
                % If step not changed too much keep Jacobian and reuse LU
                SkipLU = ( ( Theta <= OPTIONS.ThetaMin ) && (Hratio >= OPTIONS.Qmin) ...
                    && (Hratio <= OPTIONS.Qmax) );
                if ( ~SkipLU )
                    H = Hnew;
                end
            end
            SkipJac = false;
            
            % for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Accepted step. Time = ', num2str(T), ';    Stepsize = ', num2str(H) ];
                disp(str);
            end;

        else % reject  
            if ( FirstStep || Reject )
                H = OPTIONS.FacRej*H;
            else
                H = Hnew;
            end
            Reject = true;
            SkipJac = true;
            SkipLU = false;
            if ( ISTATUS.Nacc >= 1 )
                ISTATUS.Nrej = ISTATUS.Nrej + 1;
            end

            % for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Rejected step. Time = ', num2str(T), ';    Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        end % accept/reject        
    end % Tloop
    
    Ierr = 1;
    
    % Deallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        Tout(TYindex:OPTIONS.Max_no_steps) = [];
        Yout(:,TYindex:OPTIONS.Max_no_steps) = []; 
    else
        Tout = T;
        Yout = Y;
    end
    Tout = Tout';
    Yout = Yout';
            
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