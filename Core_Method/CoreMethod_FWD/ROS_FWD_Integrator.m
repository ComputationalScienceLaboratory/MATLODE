%% ROS_FWD_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ROS_FWD_Integrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )
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
% Rosenbrock forward core method.
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
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ROS_FWD_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )

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
    Roundoff = eps/2;
    
    DeltaMin = 1e-5;
    
   if ( OPTIONS.storeCheckpoint == true )
        Yout = zeros(NVAR,OPTIONS.Max_no_steps);
        Tout = zeros(OPTIONS.Max_no_steps,1);
        TYindex = 1;
        Yout(:,TYindex) = Y;
        Tout(TYindex,1) = Tinitial; 
    end    
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ros_A ros_C ros_M ros_E
    global ros_Alpha ros_Gamma
    global ros_NewF
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');    
    
    K = zeros(NVAR*Coefficient.NStage,1);
    L = zeros(OPTIONS.NADJ*Coefficient.NStage,1);
    
    quadrature = OPTIONS.Quadrature();    

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
        
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H );
        end
        if ( ( ( T+0.1*H ) == T) || ( H <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H );
        end
        
        % Limit H if necessary to avoid going beyond Tfinal
        H = min( H, abs( Tfinal-T ) );
        
        % Compute the function at current time
        Fcn0 = OdeFunction( T, Y );
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        
        % Compute the function derivative with respect to T
        dFdT=0;
        if ( ~OPTIONS.Autonomous )
            [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn0, OdeFunction, ISTATUS );
        end
        
        % Compute the Jacobian at current time
        if ( ~OPTIONS.MatrixFree )
            fjac = OPTIONS.Jacobian(T,Y);
            ISTATUS.Njac = ISTATUS.Njac + 1;
        else
            if( ~isempty(OPTIONS.Jacobian) )
                if( nargin(OPTIONS.Jacobian) == 3 )
                    fjac = @(vee)OPTIONS.Jacobian(T,Y,vee);
                elseif( nargin( OPTIONS.Jacobian )== 2 )
                    Jac = OPTIONS.Jacobian(T,Y);
                    %ISTATUS.Njac = ISTATUS.Njac + 1;
                    fjac = @(vee)(Jac*vee);
                else
                    error('Jacobian function takes a fucked up number of variables.')
                end
            else
               % Fcn0 = OdeFunction(T,Y);
	       % ISTATUS.Nfun= ISTATUS.Nfun+1;
                normy = norm(Y);
                fjac = @(v)Mat_Free_Jac(T,Y,v,OdeFunction,Fcn0,normy);
            end
        end
        
        % Repeat step calculation until current step accepted
        accepted = false;
        gmresFlag = 0;
        singCount = 0;
        while ( ~accepted ) % accepted
            
            if ( ~OPTIONS.MatrixFree )
                [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( NVAR, H, Direction, ros_Gamma(1), fjac, ISTATUS );
            else
                hgaminv = 1.0/(Direction*H*ros_Gamma(1));
                e = @(v)(hgaminv*v - fjac(v));
                ISING = 0;
            end
            
            if ( ISING ~= 0 || singCount >= 15 ) % More than 5 consecutive failed decompositions
                error('Matrix is repeatedly singular');
            end
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for istage=1:Coefficient.NStage % stages
                
                % Current istage offset. Current istage vector is K(ioffset+1:ioffset_NVAR)
                ioffset = NVAR*(istage-1);
                
                % For the 1st stage the function has been computed previously
                if ( istage == 1 )
                    Fcn = Fcn0;
                    if ( adjStackFlag )
                        Ystage = Y;
                        Ynew = Y;
                    end
                elseif ( ros_NewF(istage) )
                    Ynew = Y;
                    for j=1:istage-1;
                        Ynew = Ynew + ros_A((istage-1)*(istage-2)/2+j)*K(NVAR*(j-1)+1:NVAR*(j-1)+NVAR);
                    end
                    Tau = T + ros_Alpha(istage)*Direction*H;
                    Fcn = OdeFunction( Tau, Ynew );
                    ISTATUS.Nfun = ISTATUS.Nfun + 1;
                end

                if ( adjStackFlag )
                    if ( istage > 1 )
                        Ystage(ioffset+1:ioffset+NVAR) = Ynew;
                    end
                end

                K(ioffset+1:ioffset+NVAR) = Fcn;
                for j=1:istage-1
                    HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H);
                    K(ioffset+1:ioffset+NVAR) = K(ioffset+1:ioffset+NVAR) + HC*K(NVAR*(j-1)+1:NVAR*(j-1)+NVAR);
                end
                if ( ( ~OPTIONS.Autonomous ) && ( ros_Gamma(istage) ~= 0.0 ) )
                    HG = Direction*H*ros_Gamma(istage);
                    K(ioffset+1:ioffset+NVAR) = K(ioffset+1:ioffset+NVAR) + HG*dFdT;
                end
                if ( ~OPTIONS.MatrixFree )                
			if ( OPTIONS.LU )
			    %%%%% LU Decomp NEW %%%%%
			    [L,U,P] = lu(e);
			    K(ioffset+1:ioffset+NVAR) = U\(L\(P*K(ioffset+1:ioffset+NVAR)));
			    %%%%%%%%%%%%%%%%%%%%%%%%%
			else 
			    % Solve the system
			    K(ioffset+1:ioffset+NVAR) = e\K(ioffset+1:ioffset+NVAR);
			end                             
                else
                    [ tempK, gmresFlag, ~, iter] = ...
                        gmres(e, K(ioffset+1:ioffset+NVAR), ...
                        OPTIONS.GMRES_Restart,...
                        OPTIONS.GMRES_TOL,min(OPTIONS.GMRES_MaxIt,NVAR), OPTIONS.GMRES_P);
                    
                    if ( ~isempty(OPTIONS.GMRES_Restart) )
                         vecCount = iter(2) + (OPTIONS.GMRES_Restart - 1)*iter(1);
                    else
                         vecCount = iter(2);
                    end
                    ISTATUS.Njac =  ISTATUS.Njac + vecCount;
                    
                    if( gmresFlag ~= 0 )
                        resvec = abs(e(tempK) - K(ioffset+1:ioffset+NVAR));
                        scalar = OPTIONS.AbsTol + OPTIONS.RelTol.*abs(K(ioffset+1:ioffset+NVAR));
                        if (norm(resvec./scalar) > sqrt(NVAR))
                            singCount = singCount + 1;
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
                        else
                            gmresFlag = 0;
                        end
                    end
                    K(ioffset+1:ioffset+NVAR) = tempK;
                end

	        ISTATUS.Nsol = ISTATUS.Nsol + 1;
                
            end % stages
            
            if( gmresFlag ~= 0 )
                H = 0.5*H;
                continue;
            end
            
            Ynew = Y;
            for j=1:Coefficient.NStage
                Ynew = Ynew + ros_M(j)*K(NVAR*(j-1)+1:NVAR*(j-1)+NVAR);
            end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error Estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute the error estimation
            Yerr = zeros(NVAR,1);
            for j=1:Coefficient.NStage
                Yerr = Yerr + ros_E(j)*K(NVAR*(j-1)+1:NVAR*(j-1)+NVAR);
            end

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
                
                if ( adjStackFlag )
                    stack_ptr = rosPush(NVAR,T,H,Ystage,K,OPTIONS,Coefficient,stack_ptr);
                end
                
                if ( adjQuadFlag )
                    RY = OPTIONS.DRDY(T,Y);
                    R = OPTIONS.QFun(T,Y);
                    if ( ~OPTIONS.Autonomous )
                        Delta = sqrt(Roundoff)*max(DeltaMin,abs(T));
                        dRdT = OPTIONS.QFun(T+Delta,Y);
                        dRdT = dRdT - R;
                        dRdT = dRdT*(1/Delta);
                    end
                    
                    % Compute the stages for the quadrature term
                    for istage=1:Coefficient.NStage
                        ioffset = OPTIONS.NADJ*(istage-1);
                        if ( istage == 1 )
                            % do nothing
                        elseif ( ros_NewF(istage) )
                            Tau = T + ros_Alpha(istage)*Direction*H;
                            R = OPTIONS.QFun(Tau,Ystage(NVAR*(istage-1)+1:NVAR*(istage-1)+1+NVAR-1));
                        end
                        
                        L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) = R;
                        for j=1:istage-1
                            HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H);
                            L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) = L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) + HC*L(OPTIONS.NADJ*(j-1)+1:(OPTIONS.NADJ*(j-1)+1)+(OPTIONS.NADJ-1));
                        end
                        
                        if ( ~OPTIONS.Autonomous && (ros_Gamma(istage) ~= 0) )
                            HG = Direction*H*ros_Gamma(istage);
                            L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) = L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) + HG*dRdT;
                        end
                        
                        for i=1:OPTIONS.NADJ
                            for j=1:NVAR
                                L(ioffset+i) = L(ioffset+i)+RY(j,i)*K(NVAR*(istage-1)+j);
                            end
                        end
                        L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1)) = L(ioffset+1:(ioffset+1)+(OPTIONS.NADJ-1))*Direction*H*ros_Gamma(1);
                    end
                    
                % Compute the new quadrature solution
                    for j=1:Coefficient.NStage
                        quadrature = quadrature + ros_M(j)*L(OPTIONS.NADJ*(j-1)+1:(OPTIONS.NADJ*(j-1)+1)+(OPTIONS.NADJ-1));
                    end                    
                end
                
                % Update time
                T = T + Direction*H;          
                
                % Update solution
                Y = Ynew;
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;                
                    Yout(:,TYindex) = Y;
                    TYindex = TYindex + 1;                
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
        Tout(TYindex:OPTIONS.Max_no_steps) = [];
        Yout(:,TYindex:OPTIONS.Max_no_steps) = [];       
    else
        Tout = T;
        Yout = Y;
    end
    Yout = transpose(Yout);
    
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