%% ROK_FWD_Integrator
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
% [1] Tony D'Augustine and Adrian Sandu. MATLODE.
%
% [2] P. Tranquilli and A. Sandu, 'Rosenbrock-Krylov Methods for Large
%     Systems of Differential Equations,' SIAM Journal on Scientific 
%     Computing 36(3), A1313–A1338 
%
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ROK_FWD_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag )

    % Force initial value matrix to be 1 X N.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = Y';
    end  

    % Get Problem Size
    NVAR = max(size(Y));

    % Initialize time variables
    Tinitial = Tspan(1);
    Tfinal = Tspan(2);
        
    % Roundoff error specific to machine
    Roundoff = eps/2;
    
    DeltaMin = 1e-5;
    
    if ( OPTIONS.storeCheckpoint == true )
        [ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );     
%         Yout = zeros(NVAR,OPTIONS.Max_no_steps);
%         Tout = zeros(OPTIONS.Max_no_steps,1);
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
    global alpha gamma c gam ROK_E b
    global istage    % Added by Arash to test JacVec convergence
    istage=1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');    
    
    stack_ptr = 0;
    
    
    K = zeros(NVAR,Coefficient.NStage);
%     lambda = zeros(OPTIONS.NBasisVectors,Coefficient.NStage);
    
    
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
        if ( ~OPTIONS.MatrixFree && ~isempty(OPTIONS.Jacobian) )
            fjac = OPTIONS.Jacobian(T,Y);
            ISTATUS.Njac = ISTATUS.Njac + 1;
        else
            if( ~isempty(OPTIONS.Jacobian) )
                fjac = @(vee)OPTIONS.Jacobian(T,Y,vee);
            else
                normy = norm(Y);
                fjac = @(vee)Mat_Free_Jac(T,Y,vee,OdeFunction,Fcn0,normy);
            end
        end
        
        % Repeat step calculation until current step accepted
        accepted = false;
        enrich = false;
        while ( ~accepted ) % accepted
            
            if(~enrich)
                [Varn, Harn, w, HEnrich] = Arnoldi_MF(fjac, Fcn0, dFdT, NVAR, OPTIONS.NBasisVectors, OPTIONS.MatrixFree); % M should be an option!!!!!!!1
                M = OPTIONS.NBasisVectors;
            else
%                 [Varn, Harn, w, HEnrich] = Arnoldi_MF_Enrich(NVAR, M, Varn, Harn, w, HEnrich, Yerr);
%                 M = M + 1; % The increment needs to be the number of vectors added.  Here it is 1, since Yerr is (N,1);
            end
            
            lambda = zeros(M,Coefficient.NStage);
            
            [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS );
            
            if ( ISING ~= 0 ) % More than 5 consecutive failed decompositions
                keyboard
                Ierr = fatOde_ROS_ErrorMessage( -8, T, H );
                return;
            end
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for istage=1:Coefficient.NStage % stages
                
                insum = Y;
                outsum = zeros(M, 1);
                if (istage > 1)
                    for j = 1:istage-1
                        insum = insum + alpha(istage,j)*K(:,j);
                        outsum = outsum + gamma(istage,j)*Harn*lambda(:,j);
                    end
                    locF = OdeFunction(T + H*c(istage), insum);
                else
                    locF = Fcn0;
                end
                phi = transpose(Varn)*locF + w;
                locLHS = eye(M)-H*gam*Harn; 
                locRHS = H*(phi + outsum);
                lambda(:,istage) = locLHS\locRHS;
                K(:,istage) = Varn*lambda(:,istage) + H*(locF - Varn*phi);
                
            end % stages
            Ynew = Y;
            Yerr = zeros(NVAR,1);
            
            for j=1:Coefficient.NStage
               Yerr = Yerr + ROK_E(j)*K(:,j);
               Ynew = Ynew + b(j)*K(:,j);
            end
           % display(sprintf('K(20791,:)= %f %f %f %f %f %f %f', K(20791,:)));          

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
                enrich = false;
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
                enrich = true;
                
            end
        end
    end
    
    % Successful return
    Ierr = 1;
    
    % Decallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
%         Tout(TYindex:OPTIONS.Max_no_steps) = [];
%         Yout(:,TYindex:OPTIONS.Max_no_steps) = [];      
        Tout(ISTATUS.Nacc+1:end,:) = [];
        Yout(:,ISTATUS.Nacc+1:end) = []; 
    else
        Tout = T;
        Yout = Y;
    end
    Yout = transpose(Yout);
    
return;

function [Varn, Harn, w, HEnrich] = Arnoldi_MF(J, f, dFdT, N, M, matrixFree)

Varn = zeros(N, M);
Harn = zeros(M, M);
w = zeros(M, 1);

w(1) = 1;
beta = sqrt(f'*f + w(1)^2);
Varn(:, 1) = f/beta;
w(1) = w(1)/beta;


for i = 1:M
    if ( matrixFree )
        zeta = J(Varn(:, i)) + dFdT*w(i);
    else 
        zeta = J*Varn(:, i) + dFdT*w(i);
    end
    
    xi = 0;

    for j = 1:i
        Harn(j,i) = zeta'*Varn(:,j) + xi*w(j);
        zeta = zeta - Harn(j,i)*Varn(:,j);
        xi = xi - Harn(j,i)*w(j);
    end
%     bignorm = sqrt(zeta'*zeta + xi^2);
%     if bignorm/tau <= .25
%         for j = 1:i
%             rho = zeta'*Varn(:,j) + xi*w(j);
%             zeta = zeta - rho*Varn(:,j);
%             xi = xi - rho*w(j);
%             Harn(j,i) = Harn(j,i) - rho;
%         end
%     end
    if i < M
       bignorm = sqrt(zeta'*zeta + xi^2);
       Harn(i+1,i) = bignorm;
       Varn(:, i+1) = zeta/Harn(i+1,i);
       w(i+1) = xi/Harn(i+1,i);
    else
       bignorm = sqrt(zeta'*zeta + xi^2);
       HEnrich = bignorm;
%        Varn(:, i+1) = zeta/Harn(i+1,i);
       HEnrich = xi/HEnrich; 
    end
end

return

function [Varn, Harn, W, Henrich] = Arnoldi_MF_Enrich(N, M, ...
                                V, H, w, HEnrich, newVec)

K = size(newVec,2);
Varn = zeros(N, M+K);
Harn = zeros(M, M+K);
W = zeros(M+K, 1);

Harn(1:M,1:M) = H(:,:);
Harn(M+1,M) = HEnrich;
Varn(:,1:M) = V(:,:);
W(1:M) = w(:);

%keyboard
for i = M+1:M+K
    
    zeta = newVec(:,i-M);
    zeta = zeta/norm(zeta);
    xi = 0;
    tau = sqrt(zeta'*zeta);
    
    for j = 1:i
        Harn(j,i) = zeta'*Varn(:,j) + xi*W(j);
        zeta = zeta - Harn(j,i)*Varn(:,j);
        xi = xi - Harn(j,i)*W(j);
    end
    %bignorm = sqrt(zeta'*zeta + xi^2);
%     if bignorm/tau <= .25
%         for j = 1:i
%             rho = zeta'*Varn(:,j) + xi*W(j);
%             zeta = zeta - rho*Varn(:,j);
%             xi = xi - rho*W(j);
%             Harn(j,i) = Harn(j,i) - rho;
%         end 
%     end
    bignorm = sqrt(zeta'*zeta + xi^2);
    if(i+1 > M+K);
        Henrich = bignorm;
        Varn(:, i) = zeta/Henrich;
        W(i) = xi/Henrich;
  
    else
        Harn(i+1,i) = bignorm;
        Varn(:, i) = zeta/Harn(i+1,i);
        W(i) = xi/Harn(i+1,i);
    end
    
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
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt,edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 