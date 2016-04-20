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
    
    % TODO: cleanup disgusting hack: ==================================== %
    global io_arnoldi
    OPTIONS.IOArnoldi = io_arnoldi;
    global arnoldi_tol
    OPTIONS.AdaptiveArnoldiTol = arnoldi_tol;
    global basis_block_size
    OPTIONS.BlockSize  = basis_block_size;
    global n_recycled_vectors
    OPTIONS.NRecycledVectors = n_recycled_vectors;
    % =================================================================== %
    
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
    
    if ( OPTIONS.storeCheckpoint )
        RSTATUS.S1residual = zeros(OPTIONS.ChunkSize, 1);
        RSTATUS.stepsizes = zeros(OPTIONS.ChunkSize, 1);
        if ( OPTIONS.NBasisVectors == 0 )
            RSTATUS.basisSizes = zeros(OPTIONS.ChunkSize, 1);
        end
        if ( OPTIONS.IOArnoldi )
            RSTATUS.blockSizes = zeros(OPTIONS.ChunkSize, 1);
        end
    end
    
    stack_ptr = 0;
    
    
    K = zeros(NVAR,Coefficient.NStage);
    benrich = zeros(NVAR, OPTIONS.NRecycledVectors);
    zenrich = zeros(OPTIONS.NRecycledVectors,1);
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
    M = OPTIONS.NBasisVectors;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Direction > 0 ) && ( (T-Tfinal)+Roundoff <= 0.0 ) ...
            || ( Direction < 0 ) && ( (Tfinal-T)+Roundoff <= 0.0 ) )
        
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )

            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H );
        end
        if ( ( ( T+0.1*H ) == T) || ( abs(H) <= Roundoff ) )
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
                normy = norm(Y);
                fjac = @(v)Mat_Free_Jac(T,Y,v,OdeFunction,Fcn0,normy);
            end
        end
        
%%%%%%%%%% Prepare For Vector Recycling. %%%%%%%%%%%
        if( OPTIONS.NRecycledVectors && T ~= Tspan(1) )
                R = min(OPTIONS.NRecycledVectors, M-1);
                benrich(:,1:R) = Varn(:,end-R:end-1);
                zenrich(1:R) = w(end-R:end-1);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            
        [H, Varn, Harn, hmp1, vmp1, Lh, Uh, ph, w, VtV, Lv, Uv, pv, lambda1, phi, M, residual, io_arn_flag, ISTATUS] ...
             = ROK_ArnoldiAdapt(fjac, Fcn0, dFdT, NVAR, H, Direction, gam, ...
                                OPTIONS, ISTATUS);

        if( OPTIONS.NRecycledVectors && T ~= Tspan(1) )
            [Varn, Harn, w, Lv, Uv, pv, M, ISTATUS] ...
                = ROK_ArnoldiEnrich(fjac, dFdT, Varn, w, Harn, ...
                VtV, benrich, zenrich, R, M, NVAR, io_arn_flag, ...
                OPTIONS, ISTATUS);
            
            [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                = computeFirstStage(Fcn0, Harn, Varn, w, Lv, Uv, pv, H, M, ... 
                Direction, gam, io_arn_flag, ISTATUS, OPTIONS);
            
        end
    
        lambda = zeros(M,Coefficient.NStage);
        lambda(:,1) = lambda1;

        % Repeat step calculation until current step accepted
        accepted = false;
        while ( ~accepted ) % accepted
            
%            [Varn, Harn, w, hmp1, VtV, io_arn_flag] = Arnoldi_MF(fjac, Fcn0, dFdT, NVAR, M, OPTIONS.BlockSize, OPTIONS.MatrixFree, OPTIONS.IOArnoldi);
%            if ( OPTIONS.MatrixFree )
%                ISTATUS.Njac = ISTATUS.Njac + M;
%            end
            
            
%             [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS );
%            [H, ISING, L, U, p, ISTATUS] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS );
 
%            if ( io_arn_flag )
%                [Lv, Uv, pv] = lu(VtV);
%            end
            
            % Stage 1
%            phi = transpose(Varn)*Fcn0;
%            if ( io_arn_flag )
%                phi = Uv\(Lv\phi(pv));
%            end
%            phi = phi + w;
%            locRHS = H * phi;
%            lambda(:,1) = U\(L\(locRHS(p)));
            % Compute residual.
%            residual = abs(H * gam * hmp1) * lambda(M,1)

            % DEBUG: true residual.
%            debugK = Varn * lambda(:,1) + H*(Fcn0 - Varn*phi);
%            if ( OPTIONS.MatrixFree )
%                residual_true = norm((debugK - H*gam*fjac(debugK)) - H*Fcn0)
%            else
%                residual_true = norm((speye(NVAR) - H*gam*fjac)*debugK - H*Fcn0)
%            end

            % Solution of first stage.
            if ( RejectLastH )
                [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                    = computeFirstStage(Fcn0, Harn, Varn, w, Lv, Uv, pv, H, M, ... 
                    Direction, gam, io_arn_flag, ISTATUS, OPTIONS);
                lambda(:,1) = lambda1;
            end

            K(:,1) = Varn * lambda(:,1) + H*(Fcn0 - Varn*phi);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Stage 2 - NStage
            for istage=2:Coefficient.NStage % stages
                
                insum = Y;
                outsum = zeros(M, 1);
                if (istage > 1)
                    for j = 1:istage-1
                        insum = insum + alpha(istage,j)*K(:,j);
                        outsum = outsum + gamma(istage,j)*Harn*lambda(:,j);
                    end
                    locF = OdeFunction(T + H*c(istage), insum);
                    ISTATUS.Nfun = ISTATUS.Nfun + 1;
                else
                    locF = Fcn0;
                end
                phi = transpose(Varn)*locF;
                if ( io_arn_flag )
                    phi = Uv\(Lv\phi(pv));
                end
                phi = phi + w;
%                locLHS = eye(M)-H*gam*Harn; 
                locRHS = H*(phi + outsum);
%                lambda(:,istage) = locLHS\locRHS;
                lambda(:,istage) = Uh\(Lh\(locRHS(ph)));
                K(:,istage) = Varn*lambda(:,istage) + H*(locF - Varn*phi);

                % Stage 2 residual
%                 if ( istage == 2 )
%                     externalF = locF - Varn*phi;
%                     tempLambda = - gamma(2,1) * lambda(:,1) - gam * lambda(:,2);
%                     if ( OPTIONS.MatrixFree )
%                          residual2 = abs(H) * norm(-((H * gam) * fjac(externalF)) + ((hmp1 * tempLambda(end)) * vmp1))
%                     else
%                          residual2 = abs(H) * norm(-(H * gam * (fjac * externalF)) + (hmp1 * tempLambda(end) * vmp1))
%                     end
% 
%                     if ( OPTIONS.MatrixFree )
%                         residual2_true = norm((K(:,2) - H * gam * fjac(K(:,2))) - (H * locF) - ((H * gamma(2,1)) * fjac(K(:,1))))
%                     else
%                         residual2_true = norm(((speye(NVAR) - H * gam * fjac) * K(:,2)) - (H * locF) - ((H * gamma(2,1)) * (fjac * K(:,1))))
%                     end
%                 end
                
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
                
                % Update Memory Allocation
                if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
                    [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
                    RSTATUS.S1residual = [RSTATUS.S1residual; zeros(OPTIONS.ChunkSize, 1)];
                    RSTATUS.stepsizes = [RSTATUS.stepsizes; zeros(OPTIONS.ChunkSize, 1)];
                    if ( OPTIONS.NBasisVectors == 0 )
                        RSTATUS.basisSizes = [RSTATUS.basisSizes; zeros(OPTIONS.ChunkSize, 1)];
                    end
                    if ( OPTIONS.IOArnoldi )
                        RSTATUS.blockSizes = [RSTATUS.blockSizes; zeros(OPTIONS.ChunkSize, 1)];
                    end
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;                
                    Yout(:,TYindex) = Y;
                    RSTATUS.S1residual(TYindex) = residual;
                    RSTATUS.stepsizes(TYindex) = H;
                    if ( OPTIONS.NBasisVectors == 0 )
                        RSTATUS.basisSizes(TYindex) = M;
                    end
                    if ( OPTIONS.IOArnoldi )
                        RSTATUS.blockSizes(TYindex) = OPTIONS.BlockSize;
                    end
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
                enrich = false;
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    if ( OPTIONS.NBasisVectors == 0 )
                        str = ['Accepted step. ', 'Residual = ', num2str(residual), '; Basis Size = ', num2str(M), '; Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    else
                        str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    end
                    disp(str); 
                end
                
            else % Reject step
                if ( RejectMoreH )
                    Hnew = H*OPTIONS.FacRej;
                end
                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ) );
                
                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H = Hnew;
                
                if ( ISTATUS.Nacc >= 1 ) 
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    if ( OPTIONS.NBasisVectors == 0 )
                        str = ['Rejected step. ', 'Residual = ', num2str(residual), '; Basis Size = ', num2str(M), '; Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    else
                        str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    end
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
        RSTATUS.S1residual = RSTATUS.S1residual(1:ISTATUS.Nacc);
        RSTATUS.stepsizes = RSTATUS.stepsizes(1:ISTATUS.Nacc);
        if ( OPTIONS.NBasisVectors == 0 )
            RSTATUS.basisSizes = RSTATUS.basisSizes(1:ISTATUS.Nacc);
        end
        if ( OPTIONS.IOArnoldi )
            RSTATUS.blockSizes = RSTATUS.blockSizes(1:ISTATUS.Nacc);
        end
    else
        Tout = T;
        Yout = Y;
    end
    Yout = transpose(Yout);

return;

function [H, Varn, Harn, hmp1, vmp1, Lh, Uh, ph, w, VtV, Lv, Uv, pv, lambda1, phi, M, residual, io_arn_flag, ISTATUS] ...
                = ROK_ArnoldiAdapt(J, f, dFdT, N, H, Direction, gam, OPTIONS, ISTATUS)

Lv = []; Uv = []; pv = [];
% TODO: Some logic to choose when to do incomplete orthogonalization?
io_arn_flag = OPTIONS.IOArnoldi;

% Check whether to perform adaptive selection.
M = OPTIONS.NBasisVectors;
adaptive = false;
if ( M == 0 )
    M = 4;
    adaptive = true;
end

basisMax = min(N, 200);

testIndex = [4,6,8,11,15,20,27,36,46,57,70,85,100];
Varn = zeros(N, M);
Harn = zeros(M, M);
VtV = eye(M, M);
w = zeros(M, 1);
hmp1 = 0.0;
vmp1 = zeros(N, 1);

w(1) = 1;
beta = sqrt(f'*f + w(1)^2);
Varn(:, 1) = f/beta;
w(1) = w(1)/beta;

for i = 1:basisMax
    if ( OPTIONS.MatrixFree )
        zeta = J(Varn(:, i)) + dFdT*w(i);
        ISTATUS.Njac = ISTATUS.Njac + 1;
    else 
        zeta = J*Varn(:, i) + dFdT*w(i);
    end
    
    xi = 0;
    tau = sqrt(zeta'*zeta);

    if ( io_arn_flag )
%        error('Incomplete orthogonalization is unimplemented.')
        leftside = max(i-OPTIONS.BlockSize, 0);
        orthoVectors = (leftside+1):i;
        remainingVectors = 1:leftside;
    else
        orthoVectors = 1:i;
    end

    for j = orthoVectors
        Harn(j,i) = zeta'*Varn(:,j) + xi*w(j);
        zeta = zeta - Harn(j,i)*Varn(:,j);
        xi = xi - Harn(j,i)*w(j);
    end
    bignorm = sqrt(zeta'*zeta + xi^2);
    if bignorm/tau <= .25
        for j = orthoVectors
            rho = zeta'*Varn(:,j) + xi*w(j);
            zeta = zeta - rho*Varn(:,j);
            xi = xi - rho*w(j);
            Harn(j,i) = Harn(j,i) - rho;
        end
    end
    bignorm = sqrt(zeta'*zeta + xi^2);
    hmp1 = bignorm;

    % Check residual.
    if ( (~adaptive && i == M) || (adaptive && (i > 100 ||  min(abs(testIndex - i)) == 0)) || (i == basisMax) )
%         [H, Lh, Uh, ph, ISTATUS] = ROK_PrepareMatrix( i, H, Direction, gam, Harn(1:i, 1:i), ISTATUS );
% 
%         if ( io_arn_flag )
%             [Lv, Uv, pv] = lu(VtV(1:i, 1:i), 'vector');
%         end
% 
%         % Stage 1
%         phi = transpose(Varn(:,1:i))*f;
%         if ( io_arn_flag )
%             phi = Uv\(Lv\phi(pv));
%         end
%         phi = phi + w;
%         locRHS = H * phi;
%         lambda1 = Uh\(Lh\(locRHS(ph)));
        % Compute residual.
        if ( io_arn_flag )
            [Lv, Uv, pv, sing] = lu_decomposition(VtV);
            %v_cond = cond(VtV)
            %u_diag = diag(Uv);
            %u_diag = u_diag(max(1,i-3):end)
            if ( sing )
                warning('VtV is singlular or ill-conditioned.')
            end
        end

        [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                = computeFirstStage(f, Harn(1:i,1:i), Varn(:,1:i), w, Lv, Uv, pv, H, i, Direction, gam, io_arn_flag, ISTATUS, OPTIONS);
        
        residual = abs(H * gam * hmp1 * lambda1(i));

%        fprintf('%d: %e\n', i, residual);

        % DEBUG: true residual.
%        debugK = Varn * lambda1 + H*(f - Varn*phi);
%        if ( OPTIONS.MatrixFree )
%            residual_true = norm((debugK - H*gam*J(debugK)) - H*f);
%        else
%            residual_true = norm((speye(N) - H*gam*J)*debugK - H*f);
%        end
%        fprintf('%d: %e, %e\n', i, residual, residual_true);


        if ( ~adaptive || (i == basisMax) || residual < OPTIONS.AdaptiveArnoldiTol )
            M = i;
            break;
        end
    end

    Harn(i+1,i) = bignorm;
    Varn(:, i+1) = zeta/Harn(i+1,i);
    w(i+1) = xi/Harn(i+1,i);
    if ( io_arn_flag )
        for j = orthoVectors
            VtV(i+1,j) = 0.0;
            VtV(j,i+1) = 0.0;
        end
        for j = remainingVectors
            vtv = Varn(:,i+1)'*Varn(:,j) + w(i+1)*w(j);
            VtV(i+1,j) = vtv;
            VtV(j,i+1) = vtv;
        end
        VtV(i+1,i+1) = 1.0;
    end
end
vmp1 = zeta/hmp1;
Varn = Varn(1:N, 1:M);
Harn = Harn(1:M, 1:M);

return

function [Varn, Harn, w, Lv, Uv, pv, M, ISTATUS] ...
                = ROK_ArnoldiEnrich(J, dFdT, Varn, w, Harn, VtV, b, z, r, M, N, io_arn_flag, OPTIONS, ISTATUS)

Lv = []; Uv = []; pv = [];
for i = M:M+r-1

    zeta = b(:,i-M+1);
    xi = z(i-M+1);
    
    tau = sqrt(zeta'*zeta);

    if ( io_arn_flag )
%        error('Incomplete orthogonalization is unimplemented.')
        leftside = max(i-OPTIONS.BlockSize, 0);
        orthoVectors = (leftside+1):i;
        remainingVectors = 1:leftside;
    else
        orthoVectors = 1:i;
    end

    for j = orthoVectors
        Horth = zeta'*Varn(:,j) + xi*w(j);
        zeta = zeta - Horth*Varn(:,j);
        xi = xi - Horth*w(j);
    end
    bignorm = sqrt(zeta'*zeta + xi^2);
    if bignorm/tau <= .25
        for j = orthoVectors
            rho = zeta'*Varn(:,j) + xi*w(j);
            zeta = zeta - rho*Varn(:,j);
            xi = xi - rho*w(j);
        end
    end
    bignorm = sqrt(zeta'*zeta + xi^2);

    Varn(:, i+1) = zeta/bignorm;
    w(i+1) = xi/bignorm;
    
    if ( io_arn_flag )
        for j = orthoVectors
            VtV(i+1,j) = 0.0;
            VtV(j,i+1) = 0.0;
        end
        for j = remainingVectors
            vtv = Varn(:,i+1)'*Varn(:,j) + w(i+1)*w(j);
            VtV(i+1,j) = vtv;
            VtV(j,i+1) = vtv;
        end
        VtV(i+1,i+1) = 1.0;
    end
    
end

if ( io_arn_flag )
    [Lv, Uv, pv, sing] = lu_decomposition(VtV);
    if ( sing )
        warning('VtV is singlular or ill-conditioned.')
    end
end

if( OPTIONS.MatrixFree )
    tempVec = zeros(M+r, r);
    for i = 1:r
        tempVec(:,i) = transpose(Varn)*(J(Varn(:,M+i)) + w(M+i)*dFdT);
    end
else
    tempVec = transpose(Varn)*(J*Varn(:,M+1:M+r) + dFdT*w(M+1:M+r)');    
end

Harn = [Harn, zeros(M,r); zeros(r, M + r)];

if ( io_arn_flag )
    Harn(:,M+1:M+r) = Uv\(Lv\tempVec(pv,:));
else
    Harn(:,M+1:M+r) = tempVec;
end
    

M = M + r;
Varn = Varn(1:N, 1:M);
Harn = Harn(1:M, 1:M);

return

function [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                = computeFirstStage(f, Harn, Varn, w, Lv, Uv, pv, H, M, Direction, gam, io_arn_flag, ISTATUS, OPTIONS)

[H, Lh, Uh, ph, ISTATUS] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS, OPTIONS );

% Stage 1
phi = transpose(Varn)*f;
if ( io_arn_flag )
    phi = Uv\(Lv\phi(pv));
end
phi = phi + w;
locRHS = H * phi;
lambda1 = Uh\(Lh\(locRHS(ph)));

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
