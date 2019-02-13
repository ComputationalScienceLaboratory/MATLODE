%% EXP_FWD_Integrator
%
% <html> Up: <a href="../../../Library/html/Library.html">Library</a> </html>
%
%% Syntax
%
%
%% Input Parameters
%
%
%% Output Parameters
%
%
%% Description
%
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE
%
% [2] Hong Zhang, Adrian Sandu. FATODE: a library for forward, adjoint and
%     tangent linear integration of ODEs, SIAM Journal on Scientific
%     Computing, 36(5), C504-C523, 2014.
%
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr] = PEXP_FWD_Integrator( OdeFunction1, OdeFunction2, ...
        Tspan, Y, OPTIONS, Coefficient)
    
    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end
    
    % Get Problem Size
    NVAR = max(size(Y));
    
    % Initialize time variables
    Tinitial = Tspan(1);
    Tfinal = Tspan(2);
    
    % note: based off DLAMCH('E') NOT MATLAB's EPS
    Roundoff = 1.11022302462515654E-016;
    
    DeltaMin = 1e-5;
    
    if ( OPTIONS.storeCheckpoint == true )
        [ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );
    else
        Yout = zeros(NVAR,1);
        Tout = 0;
    end
    TYindex = 1;
    
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   Global Variables
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   Initial Settings
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
       
    T = Tinitial;
    H = min( max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) ), abs(OPTIONS.Hmax) );
    if ( abs(H) <= 10.0*Roundoff )
        H = DeltaMin;
    end
    
    if ( Tfinal >= Tinitial )
        Direction = 1;
    else
        Direction = -1;
    end
    H = Direction*H;
    
    RejectLastH = false;
    RejectMoreH = false;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   Time loop
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Direction > 0 ) && ( (T-Tfinal)+Roundoff <= 0.0 ) ...
            || ( Direction < 0 ) && ( (Tfinal-T)+Roundoff <= 0.0 ) )
        
        if (any(isnan(Y)))
            error('Solution is NaN');
        end
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H );
        end
        if ( ( ( T+0.1*H ) == T) || ( H <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H );
        end
        
        % Limit H if necessary to avoid going beyond Tfinal
        H = min( H, abs( Tfinal-T ) );
        
        % Compute the function at current time
        Fcn1_0 = OdeFunction1( T, Y );
        Fcn2_0 = OdeFunction2( T, Y );
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        
        % Compute the function derivative with respect to T
        if ( ~OPTIONS.Autonomous )
            % finite difference derivative with respect to time.
            [ dF1dT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn1_0, OdeFunction1, ISTATUS );
            [ dF2dT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn2_0, OdeFunction2, ISTATUS );
        end
        
        % Compute the Jacobian at current time
        if ( ~OPTIONS.MatrixFree )
            % Jacobian wrt y as a matrix
            fjac1 = OPTIONS.Jacobian1(T,Y);
            fjac2 = OPTIONS.Jacobian2(T,Y);
            ISTATUS.Njac = ISTATUS.Njac + 1;
        else
            if(~isempty(OPTIONS.Jacobian1) && ~isempty(OPTIONS.Jacobian2))
                % Jacobian-vector product function where Jacobian is wrt y
                fjac1 = @(vee)OPTIONS.Jacobian1(T,Y,vee);
                fjac2 = @(vee)OPTIONS.Jacobian2(T,Y,vee);
            else
                normy = norm(Y);
                % FD approx of Jacobian-vector product where Jacobian is wrt y
                fjac1 = @(vee)Mat_Free_Jac(T,Y,vee,OdeFunction1,Fcn1_0,normy);
                fjac2 = @(vee)Mat_Free_Jac(T,Y,vee,OdeFunction2,Fcn2_0,normy);
            end
        end
        
        % Repeat step calculation until current step accepted
        accepted = false;
        while ( ~accepted ) % accepted          
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Stages
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if( ~OPTIONS.Autonomous )
                y  = [Y; T];
                f1 = [Fcn1_0; 1];
                f2 = [Fcn2_0; 1];
                rhsFun1 = @(why)([OdeFunction1(why(end), why(1:end-1)); 1]);
                rhsFun2 = @(why)([OdeFunction2(why(end), why(1:end-1)); 1]);
                
                if( ~OPTIONS.MatrixFree )
                    J1 = [fjac1, dF1dT; zeros(1,NVAR+1)];
                    J2 = [fjac2, dF2dT; zeros(1,NVAR+1)];
                else
                    % computes Jy where y(1:end -1) is state space, y(end)
                    % is time
                    J1 = @(why)([fjac1(why(1:end-1)) + dF1dT*why(end); 0]);
                    J2 = @(why)([fjac2(why(1:end-1)) + dF2dT*why(end); 0]);
                end                
            else
                y       = Y;
                f1      = Fcn1_0;
                f2      = Fcn2_0;
                rhsFun1 = @(why)OdeFunction1(T, why);
                rhsFun2 = @(why)OdeFunction2(T, why);
                J1      = fjac1;
                J2      = fjac2;
            end
            
            % Need to do something about symmetricity of the Jacobians
            if ~isempty(OPTIONS.IsJACSymm)
                symmjac = OPTIONS.IsJACSymm;
            else
                symmjac = false;
            end
            
                        
            % We duplicate this because each method will set its own minimum
            % number of basis vectors which may vary from method to method,
            % in case this option is not set. For instance K-methods need
            % 4-vectors for full convergence order
            if ~isempty(OPTIONS.MBasisVectors)
                [ynew, yerr, ISTATUS] = OPTIONS.OneStepIntegrator(y,H, rhsFun1, ...
                    rhsFun2, J1, J2, f1, f2, OPTIONS.MatrixFree, OPTIONS.NBasisVectors, ISTATUS, ...
                    OPTIONS.AbsTol, OPTIONS.RelTol, OPTIONS.Adaptive_Krylov, ...
                    symmjac, OPTIONS.NReactants, OPTIONS.Autonomous, OPTIONS.MBasisVectors);
            else
                [ynew, yerr, ISTATUS] = OPTIONS.OneStepIntegrator(y,H,rhsFun1, ...
                    rhsFun2, J1, J2, f1, f2, OPTIONS.MatrixFree, OPTIONS.NBasisVectors, ISTATUS, ...
                    OPTIONS.AbsTol, OPTIONS.RelTol, OPTIONS.Adaptive_Krylov, ...
                    symmjac, OPTIONS.NReactants, OPTIONS.Autonomous);
            end
            
            if( ~OPTIONS.Autonomous )
                Ynew = ynew(1:end-1);
                Yerr = yerr(1:end-1);
            else
                Ynew = ynew;
                Yerr = yerr;
            end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Error Estimation
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute error norm
            SCAL = OPTIONS.AbsTol + OPTIONS.RelTol.*max(abs(Y),abs(Ynew));
            Err = max(sqrt(sum((Yerr./SCAL).^2)/NVAR),1d-10);
            
            % New step size is bounded by FacMin <= Hnew/H <= FacMax
            Fac = min( OPTIONS.FacMax, max( OPTIONS.FacMin, ...
                OPTIONS.FacSafe/Err^(1.0/Coefficient.ELO) ) );
            Hnew = H*Fac;
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Accept/Reject
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Check the error magnitude and adjust the step size
            ISTATUS.Nstp = ISTATUS.Nstp + 1;
            if ( ( Err <= 1.0 ) || ( H <= OPTIONS.Hmin ) )
                ISTATUS.Nacc = ISTATUS.Nacc + 1;
                
                % Update time
                T = T + Direction*H;
                
                % Update solution
                Y = Ynew;
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;
                    Yout(:,TYindex) = Y;
                    TYindex = TYindex + 1;
                end
                
                % Update Memory Allocation
                if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
                    [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end
                
                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ) );
                if ( RejectLastH ) % No step size increase after a rejected step
                    Hnew = min( Hnew, H );
                end
                
                % Last T and H
                RSTATUS.Nhexit = H;
                RSTATUS.Nhnew = Hnew;
                RSTATUS.Ntexit = T;
                
                RejectLastH = false;
                RejectMoreH = false;
                H = Hnew;
                accepted = true;
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end
                
            else % Reject step
                if ( RejectMoreH )
                    Hnew = H*OPTIONS.FacRej;
                end
                
                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H = Hnew;
                
                if ( ISTATUS.Nacc >= 1 )
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end
                
            end
        end
    end
    
    % Successful return
    Ierr = 1;
    
    % Decallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        Tout(ISTATUS.Nacc+1:end,:) = [];
        Yout(:,ISTATUS.Nacc+1:end) = [];
    else
        Tout = T;
        Yout = Y;
    end
    Yout = Yout';
    
    % Plot norm of errors
    %    figure();
    %semilogy(Yerr_nrm(1, :),'b--o')
    % hold on
    % semilogy(Yerr_nrm(2, :),'g--x')
    % semilogy(Yerr_nrm(3, :),'r*')
    % hold off
    % title(strcat(num2str(OPTIONS.Method),':',num2str(int8(-log10(OPTIONS.AbsTol)))))
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
