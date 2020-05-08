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
%     global istage    % Added by Arash to test JacVec convergence
    istage=1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
    
    if ( OPTIONS.storeCheckpoint )
        RSTATUS.S1residual = zeros(OPTIONS.ChunkSize, 1);
        RSTATUS.stepsizes = zeros(OPTIONS.ChunkSize, 1);
    end
    if ( OPTIONS.NBasisVectors == 0 )
        RSTATUS.basisSizes = zeros(OPTIONS.ChunkSize, 1);
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
        else
            dFdT = zeros(NVAR,1);
        end
        
        % Compute the Jacobian at current time
        [ fjac, ISTATUS ] = EvaluateJacobian(T, Y, Fcn0, OdeFunction, OPTIONS, ISTATUS);
        
        % Setup adjoint vector products
        if ( OPTIONS.BiOrthogonalLanczos )
            if ( ~OPTIONS.MatrixFree )
                fjact = @(v)fjac'*v;
            else
                if ( ~isempty(OPTIONS.JacobianAdjointVec) )
                    fjact = @(v)OPTIONS.JacobianAdjointVec(T,Y,v);
                elseif ( nargin(OPTIONS.Jacobian) == 2 )
                    fjact = @(v)Jac'*v;
                else
                    error('Using biorthogonal Lanczos requires some form of Jacobian Adjoint product.');
                end
            end
        end
        
        if ( OPTIONS.BiOrthogonalLanczos )
            [Varn, Wlzs, Harn, vt, wt, M, residual, H, Lh, Uh, ph, lambda1, phi, ISTATUS] = ...
              ROK_LanczosBiorthogonalize(fjac, dFdT, fjact, Fcn0, M, NVAR, H, Direction, gam, OPTIONS, ISTATUS);
        else

%%%%%%%%%% Prepare For Vector Recycling. %%%%%%%%%%%
            if( OPTIONS.NRecycledVectors && T ~= Tspan(1) )
                R = min(OPTIONS.NRecycledVectors, M-1);
                benrich(:,1:R) = Varn(:,end-R:end-1);
                zenrich(1:R) = w(end-R:end-1);
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            
            [H, Varn, Harn, hmp1, vmp1, wmp1, Lh, Uh, ph, vt, lambda1, phi, M, residual, ISTATUS] ...
                 = ROK_ArnoldiAdapt(fjac, Fcn0, dFdT, NVAR, H, Direction, gam, ...
                                    OPTIONS, ISTATUS);
    
            if( OPTIONS.NRecycledVectors && T ~= Tspan(1) )
                [Varn, Harn, w, M, ISTATUS] ...
                    = ROK_ArnoldiEnrich(fjac, dFdT, Varn, w, Harn, ...
                    benrich, zenrich, R, M, NVAR, OPTIONS, ISTATUS);
                
                [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                    = computeFirstStage(Fcn0, Harn, Varn, [], vt, [], H, M, ... 
                    Direction, gam, ISTATUS, OPTIONS);
                
            end

            Wlzs = []; wt = [];
        
        end

        lambda = zeros(M,Coefficient.NStage);
        lambda(:,1) = lambda1;

        % Repeat step calculation until current step accepted
        accepted = false;
        while ( ~accepted ) % accepted
            
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
                    = computeFirstStage(Fcn0, Harn, Varn, Wlzs, vt, wt, H, M, ... 
                    Direction, gam, ISTATUS, OPTIONS);
                lambda(:,1) = lambda1;
            end

            K(:,1) = Varn * lambda(:,1); % + H*(Fcn0 - Varn*phi);

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

                    if ( OPTIONS.BasisExtended && OPTIONS.BiOrthogonalLanczos )
%                         fprintf('OldM = %d\n', M);
                        wnew = double(~OPTIONS.Autonomous); % no time derivative if problem is autonomous.
                        [Varn, Wlzs, Harn, vt, wt, M] = ...
                            ROK_LanczosBiorthogonalEnrichBy1(fjac, dFdT, fjact, Wlzs, wt, Varn, vt, Harn, locF, wnew, M, NVAR, OPTIONS);
%                         fprintf('NewM = %d\n', M);
%                         [H, Lh2, Uh2, ph2, ISTATUS] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS, OPTIONS );
                        hm = Harn(M,1:(M-1));
                        Lh(M,:) = -H*gam*hm/Uh;
                        Lh(:,M) = [zeros(M-1,1); 1.0];
                        ph(M) = M;
                        Uh(M,1:(M-1)) = zeros(1,M-1);
                        hm = -H*gam*Harn(1:M,M); hm(M) = 1 + hm(M);
                        Uh(1:M,M) = Lh\hm(ph);
                        lambda = [lambda; zeros(1, size(lambda,2))];
                        outsum = [outsum; 0.0];
%                         ImHGH = eye(M) - H*gam*Harn;
%                         assert(norm(Lh * Uh - ImHGH(ph,:))/norm(ImHGH) < 1e-13)
%                         assert(size(Varn,2) == M && size(Harn,1) == M && size(vt,1) == M && size(lambda,1) == M && size(outsum,1) == M)
%                         assert(size(Lh,1) == M && size(Uh,1) == M && size(ph,2) == M)
                    elseif ( OPTIONS.BasisExtended )
%                        fprintf('OldM = %d\n', M);
                        wnew = double(~OPTIONS.Autonomous); % no time derivative if problem is autonomous.
                        [Varn, Harn, vt, M, ISTATUS] ...
                            = ROK_ArnoldiEnrich(fjac, dFdT, Varn, vt, Harn, ...
                            locF, wnew, 1, M, NVAR, OPTIONS, ISTATUS);
%                        fprintf('NewM = %d\n', M);
%                        [H, Lh2, Uh2, ph2, ISTATUS] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS, OPTIONS );
                        hm = Harn(1:(M-1),M);
                        Uh(:,M) = Lh\(-H*gam*hm(ph));
                        Uh(M,:) = [zeros(1,M-1), 1.0 - H*gam*Harn(M,M)];
                        Lh = [Lh, zeros(M-1,1); zeros(1,M-1), 1.0];
                        ph(M) = M;
                        lambda = [lambda; zeros(1, size(lambda,2))];
                        outsum = [outsum; 0.0];
%                         ImHGH = eye(M) - H*gam*Harn;
%                         assert(norm(Lh * Uh - ImHGH(ph,:))/norm(ImHGH) < 1e-13)
%                         assert(size(Varn,2) == M && size(Harn,1) == M && size(vt,1) == M && size(lambda,1) == M && size(outsum,1) == M)
%                         assert(size(Lh,1) == M && size(Uh,1) == M && size(ph,2) == M)
                    end
                else
                    locF = Fcn0;
                end
                
                if ( OPTIONS.BiOrthogonalLanczos )
                    phi = transpose(Wlzs)*locF + wt;
                else
                    phi = transpose(Varn)*locF + vt;
                end
%                locLHS = eye(M)-H*gam*Harn; 
                locRHS = H*(phi + outsum);
%                lambda(:,istage) = locLHS\locRHS;
                lambda(:,istage) = Uh\(Lh\(locRHS(ph)));
                K(:,istage) = Varn*lambda(:,istage) + H*(locF - Varn*phi);

                % Stage 2 residual
%                 if ( istage == 2 )
%                     tempLambda = - gamma(2,1) * lambda(:,1) - gam * lambda(:,2);
%                     if ( OPTIONS.BasisExtended )
%                         if ( OPTIONS.MatrixFree )
%                             jvhat = fjac(Varn(:,M)) + vt(M) * dFdT;
%                         else
%                             jvhat = fjac * Varn(:,M) + vt(M) * dFdT;
%                         end
%                         res2smvec = Varn' * jvhat; 
%                         res2vec = ((-1.0 * gam * lambda(M,2)) * (jvhat - Varn*res2smvec) + ((hmp1 * tempLambda(M-1)) * vmp1));
%                         residual2 = abs(H) * sqrt(res2vec'*res2vec + ((-1.0 * gam * lambda(M,2) *(-1.0 * vt' * res2smvec)) + (hmp1 * tempLambda(M-1) * wmp1))^2)
%                     else
%                         externalF = locF - Varn*phi;
%                         if ( OPTIONS.MatrixFree )
%                              residual2 = abs(H) * norm(-((H * gam) * fjac(externalF)) + ((hmp1 * tempLambda(end)) * vmp1))
%                         else
%                              residual2 = abs(H) * norm(-(H * gam * (fjac * externalF)) + (hmp1 * tempLambda(end) * vmp1))
%                         end
%                     end
% 
%                     if ( OPTIONS.MatrixFree )
%                         residual2_true = norm((K(:,2) - H * gam * fjac(K(:,2))) - (H * locF) - ((H * gamma(2,1)) * fjac(K(:,1))))
%                     else
%                         residual2_true = norm(((speye(NVAR) - H * gam * fjac) * K(:,2)) - (H * locF) - ((H * gamma(2,1)) * (fjac * K(:,1))))
%                     end
%                 end
                
            end % stages

%   New solutions

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
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end
                if ( OPTIONS.NBasisVectors == 0 && ISTATUS.Nacc > length(RSTATUS.basisSizes) )
                    RSTATUS.basisSizes = [RSTATUS.basisSizes; zeros(OPTIONS.ChunkSize, 1)];
                end
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;                
                    Yout(:,TYindex) = Y;
                    RSTATUS.S1residual(TYindex) = residual;
                    RSTATUS.stepsizes(TYindex) = H;
                    TYindex = TYindex + 1;
                end
                if ( OPTIONS.NBasisVectors == 0 )
                    RSTATUS.basisSizes(ISTATUS.Nacc) = M;
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
                    optstr = [];
                    if ( OPTIONS.NBasisVectors == 0 )
                        optstr = [optstr, ' Residual = ', num2str(residual), '; Basis Size = ', num2str(M), ';'];
                    end
                    str = ['Accepted step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H)];
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
                    optstr = [];
                    if ( OPTIONS.NBasisVectors == 0 )
                        optstr = [optstr, ' Residual = ', num2str(residual), '; Basis Size = ', num2str(M), ';'];
                    end
                    str = ['Rejected step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H)];
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
    else
        Tout = T;
        Yout = Y;
    end
    if ( OPTIONS.NBasisVectors == 0 )
        RSTATUS.basisSizes = RSTATUS.basisSizes(1:ISTATUS.Nacc);
    end
    Yout = transpose(Yout);

return;

function [H, Varn, Harn, hmp1, vmp1, wmp1, Lh, Uh, ph, w, lambda1, phi, M, residual, ISTATUS] ...
                = ROK_ArnoldiAdapt(J, f, dFdT, N, H, Direction, gam, OPTIONS, ISTATUS)

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
w = zeros(M, 1);
hmp1 = 0.0;
vmp1 = zeros(N, 1);
residual = Inf;

if ( OPTIONS.Autonomous )
    w(1) = 0;
else
    w(1) = 1;
end
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

    orthoVectors = 1:i;

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
        % Compute residual.
        [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                = computeFirstStage(f, Harn(1:i,1:i), Varn(:,1:i), [], w, [], H, i, Direction, gam, ISTATUS, OPTIONS);

        residual = abs(H * gam * hmp1 * lambda1(i));

        %fprintf('%d: %e\n', i, residual);

        % DEBUG: true residual.
%         debugK = Varn * lambda1 + H*(f - Varn*phi);
%         if ( OPTIONS.MatrixFree )
%             residual_true = norm((debugK - H*gam*J(debugK)) - H*f);
%         else
%             residual_true = norm((speye(N) - H*gam*J)*debugK - H*f);
%         end
%         fprintf('%d: %e, %e\n', i, residual, residual_true);

        if ( ~adaptive || (i == basisMax) || residual < OPTIONS.AdaptiveArnoldiTol )
            M = i;
            break;
        end
    end

    Harn(i+1,i) = bignorm;
    Varn(:, i+1) = zeta/Harn(i+1,i);
    w(i+1) = xi/Harn(i+1,i);
end
wmp1 = xi/hmp1;
vmp1 = zeta/hmp1;
Varn = Varn(1:N, 1:M);
Harn = Harn(1:M, 1:M);

return

% =========================================================================== %

function [Varn, Harn, w, M, ISTATUS] ...
                = ROK_ArnoldiEnrich(J, dFdT, Varn, w, Harn, b, z, r, M, N, OPTIONS, ISTATUS)

for i = M:M+r-1

    zeta = b(:,i-M+1);
    xi = z(i-M+1);
    
    tau = sqrt(zeta'*zeta);

    orthoVectors = 1:i;

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
Harn(:,M+1:M+r) = tempVec;

M = M + r;
Varn = Varn(1:N, 1:M);
Harn = Harn(1:M, 1:M);

return

% =========================================================================== %

function [V, W, T, vt, wt, M, residual, H, Lt, Ut, pt, lambda1, phi, ISTATUS] = ...
    ROK_LanczosBiorthogonalize(fjac, dFdT, Jt, f, M, N, H, Direction, gam, OPTIONS, ISTATUS)

if ( OPTIONS.MatrixFree )
    J = @(v)fjac(v);
else
    J = @(v)fjac*v;
end

% Check whether to perform adaptive selection.
if ( OPTIONS.RecycleBasisSize )
    oldM = M;
else
    oldM = 4;
end
M = OPTIONS.NBasisVectors;
adaptive = false;
if ( M == 0 )
    M = 4;
    adaptive = true;
end

basisMax = min(N, 200);

testIndex = [4,6,8,11,15,20,27,36,46,57,70,85,100];
ti = 1;
V = zeros(N, M);
W = zeros(N, M);
T = zeros(M, M);
vt = zeros(M, 1);
wt = zeros(M, 1);
residual = Inf;

if ( OPTIONS.Autonomous )
  vt(1) = 0.0;
  wt(1) = 0.0;
else
  vt(1) = 1.0;
  wt(1) = 1.0;
end

tau = sqrt(f'*f + vt(1)*vt(1));
V(:,1) = (1.0/tau) * f;
W(:,1) = (1.0/tau) * f;
vt(1) = vt(1)/tau;
wt(1) = wt(1)/tau;

vhat = J(V(:,1)) + dFdT * vt(1);
x = 0.0;
what = Jt(W(:,1)); % + zeros(...) * wt(1);
if ( OPTIONS.Autonomous )
  y = 0.0;
else
  y = dFdT'*W(:,1); % + 0 * wt(1);
end

alpha = vhat'*W(:,1) + x*wt(1);

vhat = vhat - alpha * V(:,1);
x = x - alpha * vt(1);
what = what - alpha * W(:,1);
y = y - alpha * wt(1);

beta = vhat'*what + x*y;
delta = sqrt(vhat'*vhat + x*x);
%delta = sqrt(abs(beta));
beta = beta / delta;

W(:,2) = (1.0/beta) * what;
wt(2) = (1.0/beta) * y;
V(:,2) = (1.0/delta) * vhat;
vt(2) = (1.0/delta) * x;
T(1,1) = alpha;
T(1,2) = beta;
T(2,1) = delta;

for j = 2:basisMax
  vhat = J(V(:,j)) + dFdT * vt(j);
  x = 0.0;
  what = Jt(W(:,j));
  if ( OPTIONS.Autonomous )
    y = 0.0;
  else
    y = dFdT'*W(:,j);
  end

  alpha = vhat'*W(:,j) + x*wt(j);

  vhat = vhat - alpha * V(:,j) - beta * V(:,j-1);
  x = x - alpha * vt(j) - beta * vt(j-1);
  what = what - alpha * W(:,j) - delta * W(:,j-1);
  y = y - alpha * wt(j) - delta * wt(j-1);

  beta = vhat'*what + x*y;
  delta  = sqrt(vhat'*vhat + x*x);
%  delta = sqrt(abs(beta));
  beta = beta / delta;

  T(j,j) = alpha;
  
  % Check residual
%  if ( (~adaptive && j == M) || (adaptive && (j > 100 || min(abs(testIndex - j)) == 0)) || (j == basisMax) )
  if ( ~adaptive && j == M )
    [H, Lt, Ut, pt, lambda1, phi, ISTATUS] ...
                  = computeFirstStage(f, T(1:j,1:j), V(:,1:j), W(:,1:j), vt(1:j), wt(1:j), H, j, Direction, gam, ISTATUS, OPTIONS);
    residual = abs(H * gam * delta * lambda1(j));
    break;
  elseif ( (adaptive && ((j > 100 && j < basisMax) || testIndex(ti) == j)) )
   if (j < 100)
    ti = ti + 1;
   end
   if ( j > 100 || testIndex(ti) >= oldM )
    [H, Lt, Ut, pt, lambda1, phi, ISTATUS] ...
                  = computeFirstStage(f, T(1:j,1:j), V(:,1:j), W(:,1:j), vt(1:j), wt(1:j), H, j, Direction, gam, ISTATUS, OPTIONS);
    residual = abs(H * gam * delta * lambda1(j));
   end
   % DEBUG: true residual.
   %debugK = V * lambda1 + H*(f - V*phi);
   %residual_true = norm((speye(N) - H*gam*J)*debugK - H*f);
   %fprintf('%d: %e, %e\n', j, residual, residual_true);
   if ( residual < OPTIONS.AdaptiveArnoldiTol )
     M = j;
     break;
   end
  elseif ( j == basisMax )
    [H, Lt, Ut, pt, lambda1, phi, ISTATUS] ...
                  = computeFirstStage(f, T(1:j,1:j), V(:,1:j), W(:,1:j), vt(1:j), wt(1:j), H, j, Direction, gam, ISTATUS, OPTIONS);
    residual = abs(H * gam * delta * lambda1(j));
    M = j;
    break;
  end
  
  T(j,j+1) = beta;
  T(j+1,j) = delta;
  W(:,j+1) = (1.0/beta) * what;
  wt(j+1) = (1.0/beta) * y;
  V(:,j+1) = (1.0/delta) * vhat;
  vt(j+1) = (1.0/delta) * x;
  
end

V = V(1:N, 1:M);
vt = vt(1:M);
W = W(1:N, 1:M);
wt = wt(1:M);
T = T(1:M, 1:M);

return

% =========================================================================== %

function [V, W, T, vt, wt, M] = ...
    ROK_LanczosBiorthogonalEnrich(fjac, dFdT, Jt, W, wt, V, vt, T, a, at, b, bt, r, M, N, OPTIONS)

if ( OPTIONS.MatrixFree )
    J = @(v)fjac(v);
else
    J = @(v)fjac*v;
end

for i = 1:r
    
    [va, vat] = ROK_GramSchmidt(W, wt, a(:,i), at(i), M+i-1, N);
    [wb, wbt] = ROK_GramSchmidt(V, vt, b(:,i), bt(i), M+i-1, N);
    dot_vw = va'*wb + vat*wbt;
    binorm = sqrt(abs(dot_vw));
    va = sign(dot_vw)*va/binorm; vat = sign(dot_vw)*vat/binorm;
    wb = wb/binorm; wbt = wbt/binorm;
    
    jva = J(va) + dFdT * vat;
    ta = W' * jva;
    tb = V' * Jt(wb) + vt * (dFdT' * wb);
    tab = wb' * jva;
    
    V(:,M+i) = va;
    vt(M+i) = vat;
    W(:,M+i) = wb;
    wt(M+i) = wbt;
    
    T(:, M+i) = ta;
    T(M+i,1:M+i-1) = tb;
    T(M+i,M+i) = tab;
end

M = M+r;

return

function [V, W, T, vt, wt, M] = ...
    ROK_LanczosBiorthogonalEnrichBy1(fjac, dFdT, Jt, W, wt, V, vt, T, a, at, M, N, OPTIONS)

if ( OPTIONS.MatrixFree )
    J = @(v)fjac(v);
else
    J = @(v)fjac*v;
end

Wa = W'*a + wt*at;
va  = (a - V*Wa);
vat = (at - vt'*Wa);
Vaug = [V, va, a; vt', vat, at];
b = [zeros(M,1); 1; 1];
if ( OPTIONS.NBasisVectors == 0 )
    [wb, ~] = lsqr(Vaug', b, OPTIONS.AdaptiveArnoldiTol); % Use linear-system tolerance for least-squares solve.
else
    wb = Vaug'\b;
end
wbt = wb(end);
wb = wb(1:end-1);

jva = J(va) + dFdT * vat;
ta = W' * jva;
tb = V' * Jt(wb) + vt * (dFdT' * wb);
tab = wb' * jva;

V(:,M+1) = va;
vt(M+1) = vat;
W(:,M+1) = wb;
wt(M+1) = wbt;

T(1:M,M+1) = ta;
T(M+1,1:M) = tb;
T(M+1,M+1) = tab;

M = M+1;

return


% =========================================================================== %

function [ba, bat] = ROK_GramSchmidt(A, at, b, bt, M, N)

tau = sqrt(b'*b + bt*bt);

ba = b/tau;
bat = bt/tau;

for j = 1:M
    scale = ba'*A(:,j) + bat*at(j);
    ba = ba - scale*A(:,j);
    bat = bat - scale*at(j);
end
nrm = sqrt(ba'*ba + bat*bat);
if nrm <= 0.25
    for j = 1:M
        rho = ba'*A(:,j) + bat*at(j);
        ba = ba - rho*A(:,j);
        bat = bat - rho*at(j);
    end
end
ba = ba*tau;
bat = bat*tau;

return

% =========================================================================== %

function [H, Lh, Uh, ph, lambda1, phi, ISTATUS] ...
                = computeFirstStage(f, Harn, Varn, Warn, vt, wt, H, M, Direction, gam, ISTATUS, OPTIONS)

[H, Lh, Uh, ph, ISTATUS] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS, OPTIONS );

% Stage 1
if ( OPTIONS.BiOrthogonalLanczos )
    phi = transpose(Warn)*f + wt;
else
    phi = transpose(Varn)*f + vt;
end
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
