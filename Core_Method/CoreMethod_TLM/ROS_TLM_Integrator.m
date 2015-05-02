function [ Tout, Yout, Y_TLM, ISTATUS, RSTATUS, Ierr ] = ROS_TLM_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: ROS_TLM_Integrator.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name    Type
%
% Output Arguments:
%   Name        Type 
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ROS_TLM_Integrator:
%
% ROS_TLM_Integrator: INPUT ARGUMENTS
%
% ROS_TLM_Integrator: OUTPUT ARGUMENTS
%
% ROS_TLM_Integrator: GLOBAL VARIABLES
%
% ROS_TLM_Integrator: SYNTAX
%
% ROS_TLM_Integrator: EXAMPLE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ros_Alpha ros_Gamma
    global ros_A ros_C ros_M ros_E
    global ros_NewF
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local variables
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

    Ynew = zeros(NVAR,1);
    Fcn0 = zeros(NVAR,1);
    Fcn = zeros(NVAR,1);
    K = zeros(NVAR*Coefficient.NStage,1);
    Ynew_tlm = zeros(NVAR,NTLM);
    Fnc0 = zeros(NVAR,NTLM);
    Fcn_tlm = zeros(NVAR,NTLM);
    K_tlm = zeros(NVAR*Coefficient.NStage,OPTIONS.NTLM);
    dFdT = zeros(NVAR,1);
    Tmp = zeros(NVAR,1);
    Yerr = zeros(NVAR,1);
    Yerr_tlm = zeros(NVAR,1);
    Pivot = zeros(NVAR,1);
    
    % Initialize ISTATUS and RSTATUS
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
    
    % Allocate for speed
    Fcn0_tlm = zeros(NVAR,NTLM);
    Fcn_tlm = zeros(NVAR,NTLM);
    
    Roundoff = 1e-250;

    TYindex = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Requirement
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( isempty(OPTIONS.Y_TLM) ) 
        error('User defined option parameter Y_TLM is required.');
    end
    if ( isempty(OPTIONS.Hess_vec) ) 
        error('User defined option parameter Hess_vec is required.');
    end    
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Allocate Memory
    [ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );

    DeltaMin = 1d-6;

    T = Tspan(1);
    Tstart = Tspan(1);
    Tend = Tspan(2);
    RSTATUS.Nhexit = 0.0;
    H = min( max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) ), abs(OPTIONS.Hmax) );
    if ( abs(H) <= 10.0*Roundoff )
        H = DeltaMin;
    end
    
    if ( Tend >= Tstart )
        Direction = 1;
    else
        Direction = -1;
    end
    H = Direction*H;
    
    RejectLastH = false;
    RejectMoreH = false;
    Transp = false;
    
    Yout = Y;
    Tout = Tstart;
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( (Direction > 0.0) && ( (T-Tend)+Roundoff <= 0.0 ) || (Direction < 0.0 ) ...
            && ( (Tend-T)+Roundoff <= 0.0 ) )
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;   H= %f', T, H );
        end
        if ( ((T+0.1*H) == T) || (H <= Roundoff) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H );
        end
        
        % Limit H if necessary to avoid going beyond Tend
        Hexit = H;
        H = min( H, abs(Tend-T) );
        
        % Compute the function at current time
        Fcn0 = OdeFunction( T, Y );
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        
        %Compute the Jacobian at current time
        fjac = OPTIONS.Jacobian(T,Y); 
        ISTATUS.Njac = ISTATUS.Njac + 1;
                
        % Compute the TLM function at current time
        Fcn0_tlm = fjac*OPTIONS.Y_TLM;
        
        % Compute the function and Jacobian derivatives with respect to T
        if ( ~OPTIONS.Autonomous )
            [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn0, OdeFunction, ISTATUS);
            Delta = sqrt(Roundoff)*max( DeltaMin, abs(T) );
            djdt = OPTIONS.Jacobian(T+Delta,Y); % LSS_Jac2(T+Delta,Y,JAC)
            ISTATUS.Njac = ISTATUS.Njac + 1;
            djdt = djdt - fjac;
            djdt = djdt/Delta;
        end

        % Repeat step calculation until current step accepted
        accepted = false;
        while ( ~accepted )
            [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( NVAR, H, Direction, ...
                ros_Gamma(1), fjac, ISTATUS ); % EDIT 1: [ros_Gamma(:)] edited this line
            if ( ISING ) % More than 5 consecutive failed decompositions
                error('Matrix is repeatedly singular \n T= %f;  H= %f', T, H );
            end
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for istage = 1:Coefficient.NStage
                % Current istage offset. Current istage vector is K(ioffset+1:ioffset+N
                ioffset = NVAR*(istage-1);
                
                % Initialize stage solution
                Ynew = Y;
                Ynew_tlm = OPTIONS.Y_TLM;
                
                % For the 1st istage the function has been computed previously
                if ( istage == 1 )
                    Fcn = Fcn0;
                    Fcn_tlm = Fcn0_tlm; % <--- probably going to have a dimension error here
                elseif ( ros_NewF(istage) ) % istage > 1 and a new function evaluation is needed at the current istage
                    for j = 1:istage-1
                        Ynew = Ynew + ros_A((istage-1)*(istage-2)/2+j)*K(NVAR*(j-1)+1:NVAR*(j-1)+NVAR);
                        for itlm=1:NTLM                       
                            Ynew_tlm(:,itlm) = Ynew_tlm(:,itlm) + ros_A((istage-1)*(istage-2)/2+j)*K_tlm(NVAR*(j-1)+1:NVAR*(j-1)+NVAR,itlm);
                        end
                    end
                    Tau = T + ros_Alpha(istage)*Direction*H;
                    Fcn = OdeFunction(Tau,Ynew);
                    ISTATUS.Nfun = ISTATUS.Nfun + 1;
                    fjac1 = lss_jac1(Tau,Ynew,OPTIONS.Jacobian);
                    ISTATUS.Njac = ISTATUS.Njac + 1;
                    for itlm=1:NTLM
                        Fcn_tlm(:,itlm) = lss_mul_jac1(Fcn_tlm(:,itlm), Ynew_tlm(:,itlm), fjac1);
                    end
                end % if istage == 1 elseif ros_NewF(istage)
                
                K(ioffset+1:ioffset+1+NVAR-1) = Fcn;
                
                for itlm=1:NTLM
                    K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) = Fcn_tlm(:,itlm);
                end
                
                for j=1:istage-1
                    HC = ros_C( (istage-1)*(istage-2)/2+j)/(Direction*H);
                    K(ioffset+1:ioffset+1+NVAR-1) = K(ioffset+1:ioffset+1+NVAR-1) + HC*K(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1);
                    for itlm=1:NTLM
                        K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) = K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) + HC*K_tlm(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1,itlm);
                    end
                end
                
                if ( ~OPTIONS.Autonomous && (ros_Gamma(istage) ~= 0.0 ) )
                    HG = Direction*H*ros_Gamma(istage);
                    K(ioffset+1:ioffset+1+NVAR-1) = K(ioffset+1:ioffset+1+NVAR-1) + HG*dFdT;
                    for itlm=1:NTLM
                        Tmp = lss_mul_jac2(Tmp, Ynew_tlm(:,itlm), djdt ); % LSS_Mul_Jac2(Tmp,Ynew_tlm(1,itlm))
                        K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) = K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) + HG*Tmp;
                    end
                end
                    [ ISING, K(ioffset+1:ioffset+1+NVAR-1) ] = lss_solve( Transp, K(ioffset+1:ioffset+1+NVAR-1), e, ISING );
                ISTATUS.Nsol = ISTATUS.Nsol + 1;
                for itlm=1:NTLM
                    Tmp = OPTIONS.Hess_vec(T,Y,K(ioffset+1:ioffset+1+NVAR-1),OPTIONS.Y_TLM(:,itlm));
                    K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) = K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) + Tmp;
                    [ ISING, K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm) ] = lss_solve( Transp, K_tlm(ioffset+1:ioffset+1+NVAR-1,itlm), e, ISING );
                    ISTATUS.Nsol = ISTATUS.Nsol + 1;
                end
            end % stage
                     
            % Compute the new solution
            Ynew = Y;
            for j=1:Coefficient.NStage
                Ynew = Ynew + ros_M(j)*K(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1);
            end

            for itlm=1:NTLM
                Ynew_tlm(:,itlm) = OPTIONS.Y_TLM(:,itlm);
                for j=1:Coefficient.NStage
                    Ynew_tlm(:,itlm) = Ynew_tlm(:,itlm) + ros_M(j)*K_tlm(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1,itlm);
                end
            end
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Error estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Yerr = zeros(NVAR,1);
            for j=1:Coefficient.NStage
                Yerr = Yerr + ros_E(j)*K(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1);
            end

            % calculate error norm
            SCAL = ( OPTIONS.AbsTol + OPTIONS.RelTol .* max(abs(Ynew),abs(Y)) );
            Err = sum((Yerr./SCAL).^2,1);    
            Err = max(sqrt( Err/double(NVAR)), 1.0d-10);


            if ( OPTIONS.TLMTruncErr ) %<----------------original
                Yerr_tlm = zeros(NVAR,NTLM);
                for itlm=1:NTLM
                    for j=1:Coefficient.NStage
                        Yerr_tlm(:,itlm) = Yerr_tlm(:,itlm) + ros_E(j)*K_tlm(NVAR*(j-1)+1:NVAR*(j-1)+1+NVAR-1,itlm);
                    end
                end
                Err = ros_ErrorNorm_tlm( NVAR, Ynew_tlm, Yerr_tlm, Err, OPTIONS );
            end
            
            % New step size is bounded by FacMin <= Hnew/H <= FacMax
            Fac = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafe/Err^(1/Coefficient.ELO)));
            Hnew = H*Fac;
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Accept/Reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Check the error magnitude and adjust step size
            ISTATUS.Nstp = ISTATUS.Nstp + 1;
            if ( (Err <= 1.0 ) || (H <= OPTIONS.Hmin) )

                % Update Memory Allocation
                if ( ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk )
                    [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end  

                % Update time
                T = T + Direction*H;

                % Update solution
                Y = Ynew;
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,:) = T;                
                    Yout(:,TYindex) = Y;
                    TYindex = TYindex + 1;                
                end

                % Update sensitivity
                OPTIONS.Y_TLM = Ynew_tlm;
                Y_TLM = Ynew_tlm; % <<<<<<<<<<<<<<<<<<<<<< TEMP FIX

                ISTATUS.Nacc = ISTATUS.Nacc + 1;
                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ));
                if ( RejectLastH )
                    Hnew = min(Hnew,H);
                end
                RSTATUS.Nhexit = H;
                RSTATUS.Nhnew = Hnew;
                RSTATUS.Ntexit = T;
                RejectLastH = false;
                RejectMoreH = false;
                H = Hnew;
                accepted = true; % EXIT while loop
                                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);  
                end
                
            else
                if ( RejectMoreH )
                    Hnew = H*OPTIONS.FacRej;
                end
                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H = Hnew;
                if (ISTATUS.Nacc >= 1.0 )
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H) ];
                    disp(str);
                end
                
            end
        end % accepted
    end % time
    
    % Successful return
    Ierr = 1;

    % Deallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        [ Tout, Yout ] = matlOde_deallocateMemory(Tout,Yout,ISTATUS.Nacc);
    else
        Tout = T;
        Yout = Y;
    end
    
    % Temporary Fix
    Yout = transpose(Yout);    

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END: ROS_TLM_Integrator: CORE METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START: ROS_TLM_Integrator: SUBROUTINES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (START) SUBROUTINE: Rosenbrock Error Norm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Computes the "scaled norm" of the error vector Yerr
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ ros_ErrorNorm ] = ros_ErrorNorm( NVAR, Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
    
    Err = 0;
    
    for i = 1:NVAR
        Ymax = max( abs(Y(i)), abs(Ynew(i)));
        if ( VectorTol )
            Scale = AbsTol(i) + RelTol(i)*Ymax;
        else
            Scale = AbsTol(1) + RelTol(1)*Ymax;
        end
        Err = Err + ( Yerr(i)/Scale )^2;
    end
    Err = sqrt( Err/double(NVAR) );
    
    ros_ErrorNorm = max( Err, 1.0d-10 );

    % equivalent to:
%      SCAL = ( AbsTol + RelTol .* max(abs(Ynew),abs(Y)) );
%      Err = sum((Yerr./SCAL).^2,1);    
%      Err = max(sqrt( Err/double(NVAR)), 1.0d-10);
   
return; % (END) SUBROUTINE: Rosenbrock Error Norm

% (START) SUBROUTINE: Rosenbrock Error Norm Tangent Linear Model
function [ Err ] = ros_ErrorNorm_tlm( NVAR, Ynew_tlm, Yerr_tlm, ...
    Fwd_Err, OPTIONS )

%     Err = Fwd_Err;
%     for itlm = 1:OPTIONS.NTLM
%         TMP = ros_ErrorNorm( OPTIONS.NTLM, OPTIONS.Y_TLM(:,itlm), Ynew_tlm(:,itlm), Yerr_tlm(:,itlm), ...
%             OPTIONS.AbsTol_TLM(:,itlm), OPTIONS.RelTol_TLM(:,itlm), OPTIONS.ITOL );
%         Err = max( Err, TMP );
%     end
%     
%     ros_ErrorNorm_tlm = max( Err, 1.0d-10 );

    [~, tlmIndex] = max(sum(max(abs(Ynew_tlm),abs(OPTIONS.Y_TLM))));
    SCAL = OPTIONS.AbsTol_TLM(:,tlmIndex) + OPTIONS.RelTol_TLM(:,tlmIndex) .* max(abs(Ynew_tlm(:,tlmIndex)),abs(OPTIONS.Y_TLM(:,tlmIndex)));
    Err = sum((Yerr_tlm(:,tlmIndex)./SCAL).^2,1);    
    Err = max([sqrt( Err./double(NVAR)) 1.0d-10 Fwd_Err]);
        
return; % (END) SUBROUTINE: Rosenbrock Error Norm Tangent Linear Model

