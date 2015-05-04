%% RK_ADJ1_DiscreteIntegrator
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
function [ Tout, Yout, Lambda, Mu, ISTATUS, RSTATUS, Ierr ] = RK_ADJ1_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global rkA rkB rkC

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Roundoff = 1.11022302462515654E-016;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');  

    TYindex = 1;
    
    Mu = OPTIONS.Mu();    
    Lambda = OPTIONS.Lambda;
    
    WY1 = zeros(NVAR,OPTIONS.NADJ);
    WY2 = zeros(NVAR,OPTIONS.NADJ);
    WY3 = zeros(NVAR,OPTIONS.NADJ);
    
    WP1 = zeros(OPTIONS.NP,OPTIONS.NADJ);
    WP2 = zeros(OPTIONS.NP,OPTIONS.NADJ);
    WP3 = zeros(OPTIONS.NP,OPTIONS.NADJ);
    
    FPJAC1 = zeros(NVAR,OPTIONS.NP);
    FPJAC2 = zeros(NVAR,OPTIONS.NP);
    FPJAC3 = zeros(NVAR,OPTIONS.NP);
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % T1 = Tend;
    T1 = 20; %<---------------temp fix
    NewtonConverge = true;
    Reject = false;
    SkipJac = false;
    
    OPTIONS.AdjointSolve = 1; %<---- only adjoint solve available

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( stack_ptr > 0 )
        if ( ~Reject )
            % Recover checkpoints for stage values and vectors
            [ T, H, Y, Zstage, NewIt, stack_ptr ] = rkPop( NVAR, stack_ptr );
            Yout(:,TYindex) = Y;
            Tout(TYindex,1) = T;      
            
            TYindex = TYindex + 1;
            
            % Compute LU decomposition
            if ( ~OPTIONS.SaveLU )
                % Compute the Jacobian matrix
                JAC = OPTIONS.Jacobian(T,Y);
                ISTATUS.Njac = ISTATUS.Njac + 1;
                
                % Comptue the matrices E1 and E2 and their decompositions
                [ e1, e2, ISING, ISTATUS ] = RK_Decomp( NVAR, H, JAC, ISTATUS );
            end
            
            % Jacobian values at stage vectors
            TMP = Y + Zstage(1:NVAR,1);
            fjac1 = OPTIONS.Jacobian( T+rkC(1)*H, TMP);
            FPJAC1 = OPTIONS.Jacp(T+rkC(1)*H,TMP);                
            if ( adjQuadFlag )
                WY1 = OPTIONS.DRDY(T+rkC(1)*H,TMP);
                WP1 = OPTIONS.DRDP(T+rkC(1)*H,TMP);
            end
            TMP = Y + Zstage(NVAR+1:2*NVAR,1);
            fjac2 = OPTIONS.Jacobian( T+rkC(2)*H, TMP );
            FPJAC2 = OPTIONS.Jacp(T+rkC(2)*H,TMP);
            if ( adjQuadFlag )
                WY2 = OPTIONS.DRDY(T+rkC(2)*H,TMP);
                WP2 = OPTIONS.DRDP(T+rkC(2)*H,TMP);
            end
            TMP = Y + Zstage(2*NVAR+1:3*NVAR,1);
            fjac3 = OPTIONS.Jacobian( T+rkC(3)*H, TMP );
            FPJAC3 = OPTIONS.Jacp(T+rkC(3)*H,TMP);
            if ( adjQuadFlag )
                WY3 = OPTIONS.DRDY(T+rkC(3)*H,TMP);
                WP3 = OPTIONS.DRDP(T+rkC(3)*H,TMP);
            end
            ISTATUS.Njac = ISTATUS.Njac + 3;
            
            % for debugging
%             str = ['Acc: stack_ptr= ', num2str(stack_ptr),...
%                 '; Time= ', num2str(T), '; Lambda= ', num2str(Lambda(1,1)), ...
%                 ' ', num2str(Lambda(2,1)), ' ', num2str(Lambda(1,2)), ...
%                 ' ', num2str(Lambda(2,2))];
%             disp(str);
        else
            % for debugging
%             str = ['Rej: stack_ptr= ', num2str(stack_ptr),...
%                 '; Time= ', num2str(T), '; Lambda= ', num2str(Lambda(1,1)), ...
%                 ' ', num2str(Lambda(2,1)), ' ', num2str(Lambda(1,2)), ...
%                 ' ', num2str(Lambda(2,2))];
%             disp(str);
        end

% 111 CONTINUE       

% OPTIONS.AdjointSolve
%   (0,1): Solve_fixed
%       2: Solve_direct
%       3: Solve_adaptive
        if( (OPTIONS.AdjointSolve == 3 && ~NewtonConverge) || ...
                (OPTIONS.AdjointSolve == 2) )
            ax_big = lss_decomp_big( NVAR, H, fjac1, fjac2, fjac3 );
            ISTATUS.Ndec = ISTATUS.Ndec + 1;
        end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Loop for the simplified Newton iterations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for iadj=1:OPTIONS.NADJ
            U1 = zeros(NVAR,1);
            U2 = zeros(NVAR,1);
            U3 = zeros(NVAR,1);

            % Initializations for Newton iteration
            NewtonDone = false;

            % Right Hand Side - part G for all Newton iterations
            G1 = zeros(NVAR,1);
            G2 = zeros(NVAR,1);
            G3 = zeros(NVAR,1);

            TMP = transpose(fjac1)*Lambda(:,iadj);
            G1 = G1 - H*rkB(1)*TMP;
            if( adjQuadFlag )
                G1 = G1 - H*rkB(1)*WY1(:,iadj);
            end

            TMP = transpose(fjac2)*Lambda(:,iadj);
            G2 = G2 - H*rkB(2)*TMP;
            if( adjQuadFlag )
                G2 = G2 - H*rkB(2)*WY2(:,iadj);
            end

            TMP = transpose(fjac3)*Lambda(:,iadj);
            G3 = G3 - H*rkB(3)*TMP;
            if( adjQuadFlag )
                G3 = G3 - H*rkB(3)*WY3(:,iadj);
            end
            

            if( (OPTIONS.AdjointSolve == 3 && ~NewtonConverge) || ...
                    (OPTIONS.AdjointSolve == 1 || OPTIONS.AdjointSolve == 0) )

                % NewtonLoopAdj
                for NewtonIter=1:OPTIONS.NewtonMaxIt
                    
                    % Prepare the right-hand side
                    [ DU1, DU2, DU3 ] = RK_PrepareRHS_Adj(NVAR,H,Lambda(:,iadj),U1,U2,U3,G1,G2,G3,fjac1,fjac2,fjac3);
                                        
                    % Solve the linear systems
                    [ DU1, DU2, DU3, ISTATUS ] = RK_SolveTR( NVAR, H, DU1, DU2, DU3, e1, e2, ISTATUS );
                    
                    % The following code performs and adaptive number of
                    % Newton iterations for solving adjoint system
                    if( OPTIONS.AdjointSolve == 3) 
                        SCAL = fatOde_ErrorScale( NVAR, OPTIONS.ITOL, OPTIONS.AbsTol_ADJ(:,iadj), ...
                            OPTIONS.RelTol_ADJ(:,iadj), Lambda(:,iadj) );
                        
                        NewtonIncrement = sqrt( (fatOde_ErrorNorm(NVAR,SCAL,DU1)^2 + ...
                            fatOde_ErrorNorm(NVAR,SCAL,DU2)^2 + ...
                            fatOde_ErrorNorm(NVAR,SCAL,DU3)^2)/3 );

                        if ( NewtonIter == 1 ) 
                            Theta = abs(ThetaMin);
                            NewtonRate = 2.0;
                        else
                            Theta = NewtonIncrement/NewtonIncrementOld;
                            if ( Theta < 0.99 )
                                NewtonRate = Theta/(1-Theta);
                            else
                                Reject = true;
                                NewtonDone = false;
                                break;
                            end
                        end
                        NewtonIncrementOld = max(NewtonIncrement, Roundoff);
                    end

                    % Update solution
                    U1 = U1 - DU1;
                    U2 = U2 - DU2;
                    U3 = U3 - DU3;
                                                     
                    if( OPTIONS.AdjointSolve == 3 )
                        % When performing an adaptive number of iterations
                        % check the error in Newton iterations
                        NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol);
                        if ( NewtonDone && (NewtonIter > (NewIt+1)) )
                            break; % EXIT NewtonLoopAdj
                        end
                    elseif( OPTIONS.AdjointSolve == ( 0 || 1 ) )
                        if( NewtonIter > max(NewIt+1,6) )
                            break; % EXIT NewtonLoopAdj
                        end
                    end
                
                end % NewtonLoopAdj
            
            if( OPTIONS.AdjointSolve == 3 && ~NewtonDone )
                disp('Newton iterations do not converge, switching to full system.');
                NewtonConverge = false;
                Reject = true;
                % GOTO 111
            end

                % Compute Mu
                TMP = zeros(NVAR,1);
                TMP = TMP + H*rkA(1,1)*U1;
                TMP = TMP + H*rkA(2,1)*U2;
                TMP = TMP + H*rkA(3,1)*U3;
                TMP = TMP + H*rkB(1)*Lambda(:,iadj);
                V1 = FPJAC1'*TMP;
                if ( adjQuadFlag )
                    V1 = V1 + H*rkB(1)*WP1(:,iadj);
                end
                
                TMP = zeros(NVAR,1);
                TMP = TMP + H*rkA(1,2)*U1;
                TMP = TMP + H*rkA(2,2)*U2;
                TMP = TMP + H*rkA(3,2)*U3;
                TMP = TMP + H*rkB(2)*Lambda(:,iadj);
                V2 = FPJAC2'*TMP;
                if ( adjQuadFlag )
                    V2 = V2 + H*rkB(2)*WP2(:,iadj);
                end

                TMP = zeros(NVAR,1);
                TMP = TMP + H*rkA(1,3)*U1;
                TMP = TMP + H*rkA(2,3)*U2;
                TMP = TMP + H*rkA(3,3)*U3;
                TMP = TMP + H*rkB(3)*Lambda(:,iadj);
                V3 = FPJAC3'*TMP;
                if ( adjQuadFlag )
                    V3 = V3 + H*rkB(3)*WP3(:,iadj);
                end

                % Update Lambda
                Lambda(:,iadj) = Lambda(:,iadj) + U1;
                Lambda(:,iadj) = Lambda(:,iadj) + U2;
                Lambda(:,iadj) = Lambda(:,iadj) + U3;

                % Update Mu
                Mu(:,iadj) = Mu(:,iadj) + V1;
                Mu(:,iadj) = Mu(:,iadj) + V2;
                Mu(:,iadj) = Mu(:,iadj) + V3;
                            
            else % NewtonConverge == false
                X(1:NVAR) = -G1;
                X(NVAR+1:2*NVAR) = -G2;
                X(2*NVAR+1:3*NVAR) = -G3;
                
                %X = ax_big\(X'); <-delete me?
                X = ax_big\X;
                X = X';
                ISTATUS.Nsol = ISTATUS.Nsol + 1;
                
                TMP = zeros(NVAR,1);
                TMP = TMP + H*rkA(1,1)*X(1:NVAR);
                TMP = TMP + H*rkA(2,1)*X(NVAR+1:2*NVAR);
                TMP = TMP + H*rkA(3,1)*X(2*NVAR+1:3*NVAR);
                TMP = TMP + H*rkB(1)*Lambda(:,iadj);
                V1 = FPJAC1'*TMP;
                if ( adjQuadFlag )
                    V1 = V1 + H*rkB(1)*WP1(:,iadj);
                end
                
                TMP = zeros(NVAR,1);
                TMP = TMP + rkA(1,2)*X(1:NVAR);
                TMP = TMP + rkA(2,2)*X(NVAR+1:2*NVAR);
                TMP = TMP + rkA(3,2)*X(2*NVAR+1:3*NVAR);
                TMP = TMP + H*rkB(2)*Lambda(:,iadj);
                V2 = FPJAC2'*TMP;
                if ( adjQuadFlag )
                    V2 = V2 + H*rkB(2)*WP2(:,iadj);
                end
                
                TMP = zeros(NVAR,1);
                TMP = TMP + rkA(1,3)*X(1:NVAR);
                TMP = TMP + rkA(2,3)*X(NVAR+1:2*NVAR);
                TMP = TMP + rkA(3,3)*X(2*NVAR+1:3*NVAR);
                TMP = TMP + H*rkB(3)*Lambda(:,iadj);
                V3 = FPJAC3'*TMP;
                if ( adjQuadFlag )
                    V3 = V3 + H*rkB(3)*WP3(:,iadj);
                end
                
                Lambda(1:NVAR,iadj) = Lambda(1:NVAR,iadj) + X(1:NVAR) + X(NVAR+1:2*NVAR) + X(2*NVAR+1:3*NVAR);
                Mu(1:OPTIONS.NP,iadj) = Mu(1:OPTIONS.NP,iadj) + V1(1:OPTIONS.NP) + V2(1:OPTIONS.NP) + V3(1:OPTIONS.NP);
                
                if ( (OPTIONS.AdjointSolve == 3) && (iadj >= OPTIONS.NADJ) )
                    NewtonConverge = true;
                    Reject = false;
                end
                
            end
                                                         
        end
                
        T1 = T1 - H;
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
    end % end while
    
    % Successful Exit
    Ierr = 1;

return;

function [ R1, R2, R3, ISTATUS ] = RK_SolveTR( NVAR, H, R1, R2, R3, e1, e2, ISTATUS )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkAinvT rkTinv
    
    for i=1:NVAR
        x1 = R1(i)/H;
        x2 = R2(i)/H;
        x3 = R3(i)/H;
        R1(i) = rkAinvT(1,1)*x1 + rkAinvT(2,1)*x2 + rkAinvT(3,1)*x3;
        R2(i) = rkAinvT(1,2)*x1 + rkAinvT(2,2)*x2 + rkAinvT(3,2)*x3;
        R3(i) = rkAinvT(1,3)*x1 + rkAinvT(2,3)*x2 + rkAinvT(3,3)*x3;
    end
    
    R1 = (e1')\R1; 
    rhs = zeros(NVAR,1);
    for m=1:NVAR
        rhs(m) = R2(m) + R3(m)*1i;
    end
    
    % R3 = -R3;
    solution = e2'\rhs;
    R2 = real(solution);
    R3 = imag(solution);
    
    for i=1:NVAR
        x1 = R1(i);
        x2 = R2(i);
        x3 = R3(i);
        R1(i) = rkTinv(1,1)*x1 + rkTinv(2,1)*x2 + rkTinv(3,1)*x3;
        R2(i) = rkTinv(1,2)*x1 + rkTinv(2,2)*x2 + rkTinv(3,2)*x3;
        R3(i) = rkTinv(1,3)*x1 + rkTinv(2,3)*x2 + rkTinv(3,3)*x3;
    end
        
    ISTATUS.Nsol = ISTATUS.Nsol + 2;
    
return;

function [ R1, R2, R3 ] = RK_PrepareRHS_Adj(NVAR,H,Lambda,U1,U2,U3,G1,G2,G3,fjac1,fjac2,fjac3)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA

    R1 = G1;
    R2 = G2;
    R3 = G3;
    
    F1 = zeros(NVAR,1);
    F1 = F1 - H*rkA(1,1)*U1;
    F1 = F1 - H*rkA(2,1)*U2;
    F1 = F1 - H*rkA(3,1)*U3;
    TMP = fjac1'*F1;
    TMP = TMP + U1;
    R1 = R1 + TMP;
    
    F2 = zeros(NVAR,1);
    F2 = F2 - H*rkA(1,2)*U1;
    F2 = F2 - H*rkA(2,2)*U2;
    F2 = F2 - H*rkA(3,2)*U3;
    TMP = fjac2'*F2;
    TMP = TMP + U2;
    R2 = R2 + TMP;
    
    F3 = zeros(NVAR,1);
    F3 = F3 - H*rkA(1,3)*U1;
    F3 = F3 - H*rkA(2,3)*U2;
    F3 = F3 - H*rkA(3,3)*U3;
    TMP = fjac3'*F3;
    TMP = TMP + U3;
    R3 = R3 + TMP;
           
return;


% (START) SUBROUTINE: Runge Kutaa Decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the matrices E1 and E2 and their decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ e1, e2, ISING, ISTATUS ] = RK_Decomp( NVAR, H, JAC, ISTATUS )
    global rkGamma rkAlpha rkBeta

    Gamma = rkGamma/H;
    Alpha = rkAlpha/H;
    Beta = rkBeta/H;
    [ ISING, e1 ] = lss_decomp( NVAR, Gamma, JAC );
    
    % e1 = lu(e1);
    
    if ( ISING ~= 0 )
        ISTATUS.Ndec = ISTATUS.Ndec;
        e2 = e1;
        return;
    end

    [ ISING, e2 ] = lss_decomp_cmp( NVAR, Alpha, Beta, JAC );
    ISTATUS.Ndec = ISTATUS.Ndec + 1;
    
    % e2 = lu(e2);

return; % (END) SUBROUTINE: Runge Kutaa Decomposition




