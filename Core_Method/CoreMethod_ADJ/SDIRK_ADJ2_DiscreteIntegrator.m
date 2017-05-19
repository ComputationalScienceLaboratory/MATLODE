%% SDIRK_ADJ2_DiscreteIntegrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = SDIRK_ADJ2_DiscreteIntegrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, stack_ptr )
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
% |stack_ptr|: pointer for global snapshot stack
%
% |adjQuadFlag|: Adjoint quadrature stack flag
%
%% Output Parameters
% |Lambda|: Sensitivity matrix
%
% |ISTATUS|: Integer statistics 
%
% |RSTATUS|: Real statistics
%
% |Ierr|: Error flag
%
%% Description
% Singly Diagonally Implicit Runge-Kutta discrete 2 adjoint core method.
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
function [ Lambda, ISTATUS, RSTATUS, Ierr ] = SDIRK_ADJ2_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global rkA rkB rkC
    global rkGamma

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TYindex = 1;
    Lambda = OPTIONS.Lambda;
    
    SkipJac = false;
    SkupLU = false;
    Reject = false;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');    

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    while ( stack_ptr > 0 )
        % Recover checkpoints for stage values and vectors
        [ T, H, Y, Z, stack_ptr ] = sdirkPop( NVAR, Coefficient.NStage, stack_ptr );
%         Yout(:,TYindex) = Y;
%         Tout(TYindex,1) = T;        
%         TYindex = TYindex + 1;
        
        % Compute E = 1/(h*gamma)-Jac and its LU decomposition
        if ( ~OPTIONS.SaveLU )
            SkipJac = false;
            SkipLU = false;
            [ H, e1, SkipJac, SkipLU, Reject, ISTATUS, ISING ] = SDIRK_PrepareMatrix ( NVAR, H, T, Y, SkipJac, SkipLU, Reject, OPTIONS.Jacobian, ISTATUS );
            if ( ISING ~= 0 )
                warning('Matrix is repeatedly singular');
            end
        end
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
        for istage=Coefficient.NStage:-1:1
            % Jacobian of the current solution
            TMP = Y + Z(:,istage);
            fjac = OPTIONS.Jacobian(T+rkC(istage)*H, TMP);
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            if ( adjQuadFlag )
                WY = OPTIONS.DRDY(T+rkC(istage),TMP);
            end
            
            if ( OPTIONS.DirectADJ )
                HGamma = 1.0/(H*rkGamma);
                [ ISING, e1 ] = lss_decomp( NVAR, HGamma, fjac );
                ISTATUS.Ndec = ISTATUS.Ndec + 1;
                if ( ISING ~= 0 )
                    str = [ 'At stage ', num2str(istage), ' the matrix used in adjoint', ...
                        ' computation is singular' ];
                    disp(str);
                    warning('Matrix is repeatedly singular');
                end
            end
            
            for iadj=1:OPTIONS.NADJ
                % Update scaling coefficients
                SCAL = 1.0 ./ ( OPTIONS.AbsTol_ADJ + OPTIONS.RelTol_ADJ .* abs(Lambda(:,iadj)) );
                
                % Prepare the loop-independent part of the right-hand side
                % G(:) = H*Jac^T*( B(i)*Lambda + sum_j A(j,i)*Uj(:) )
                G = rkB(istage)*Lambda(:,iadj);
                if ( istage < Coefficient.NStage )
                    for j=istage+1:Coefficient.NStage
                        G = G + rkA(j,istage)*U(:,iadj,j);
                    end
                end
                
                TMP = fjac'*G;
                if ( adjQuadFlag )
                    TMP = TMP + rkB(istage)*WY(:,iadj);
                end
                G = H*TMP;
                
                if ( OPTIONS.DirectADJ )
                    trans = true;
                    [ G, ISTATUS ] = SDIRK_Solve(trans, H, NVAR, e1, G, ISTATUS );
                    U(:,iadj,istage) = G;
                    
                else
                    % Initializations for Newton iterations
                    U(:,iadj,istage) = zeros(NVAR,1,1);
                    
                    for NewtonIter=1:OPTIONS.NewtonMaxIt
                        % Prepare hte loop-dependent part of the right-hand
                        % side
                        TMP = fjac'*U(:,iadj,istage);
                        DU = U(:,iadj,istage) - (H*rkGamma)*TMP - G;
                        
                        % Solve the linear system
                        trans = true;
                        [ DU, ISTATUS ] = SDIRK_Solve(trans, H, NVAR, e1, DU, ISTATUS );
                        
                        % Check convergence of Newton iterations
                        NewtonIncrement = errorNorm(NVAR,DU,SCAL);
                        if ( NewtonIter == 1 ) 
                            Theta = abs(OPTIONS.ThetaMin);
                            NewtonRate = 2.0;
                        else
                            Theta = NewtonIncrement/NewtonIncrementOld;
                            if ( Theta < 0.99 )
                                NewtonRate = Theta/(1-Theta);
                                % Predict error at the end of Newton process
                                NewtonPredictedErr = NewtonIncrement*Theta^(OPTIONS.NewtonMaxIt-NewtonIter)/(1.0-Theta);
                                
                                % Non-convergence of Newton: predicted error too large
                                if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
                                    break; % NewtonLoop
                                end
                            else % non-convergence of Newton: Theta too large
                                break;
                            end
                        end
                        NewtonIncrementOld = NewtonIncrement;
                        
                        % Update Solution
                        U(:,iadj,istage) = U(:,iadj,istage) - DU;
                        
                        % Check error in Newton iterations
                        NewtonDone = (NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol);
                        
                        if ( (NewtonIter >= 4) && NewtonDone )
                            break; % NewtonLoop
                        end
                        
                    end
                    
                    % If Newton iterations fail emply the direct solution
                    if ( ~NewtonDone )
                        HGamma = 1.0/(H*rkGamma);
                        [ ISING, e1 ] = lss_decomp( NVAR, HGamma, fjac );
                        ISTATUS.Ndec = ISTATUS.Ndec + 1;
                        if ( ISING ~= 0 )
                            str = [ 'At stage ', num2str(istage), ' the matrix used in adjoint', ...
                                ' computation is singular' ];
                            disp(str);
                            warning('Matrix is repeatedly singular');
                        end
                        trans = true;
                        [ G, ISTATUS ] = SDIRK_Solve(trans, H, NVAR, e1, G, ISTATUS );
                        U(:,iadj,istage) = G;
                    end
                    
                end
            end  
        end
        
        % Update adjoint Solution
        for istage=1:Coefficient.NStage
            for iadj=1:OPTIONS.NADJ
                Lambda(:,iadj) = Lambda(:,iadj) + U(:,iadj,istage);
            end
        end
                
    end
    
    % Successful return
    Ierr = 1;

return;

function [ RHS, ISTATUS ] = SDIRK_Solve( Transp, H, NVAR, e1, RHS, ISTATUS )
    

    global HGammaInv rkGamma
    
    HGammaInv = 1.0/(H*rkGamma);
    RHS = RHS.*HGammaInv;
    if ( Transp == true )
        RHS = transpose(e1)\RHS;
    elseif ( Transp == false )
        RHS = e1\RHS;
    else
        disp( 'ERROR: SDIRK_Solve (Transp)' );    
    end
    ISTATUS.Nsol = ISTATUS.Nsol + 1;

return;

function [ H, e1, SkipJac, SkipLU, Reject, ISTATUS, ISING ] = SDIRK_PrepareMatrix ( NVAR, H, T, Y, SkipJac, SkipLU, Reject, Jacobian, ISTATUS )

    global HGammaInv rkGamma

    ConsecutiveSng = 0;
    ISING = 1;
    
    while ( ISING ~= 0 )
        HGammaInv = 1.0/(H*rkGamma);
        
        if ( ~SkipJac )
            fjac = Jacobian(T,Y);
            ISTATUS.Njac = ISTATUS.Njac + 1;
        end
        [ ISING, e1 ] = lss_decomp( NVAR, HGammaInv, fjac );
        ISTATUS.Ndec = ISTATUS.Ndec + 1;
        
        if ( ISING ~= 0 )
            str = [ 'MATRIX IS SINGULAR , ISING=', num2str(ISING), ...
                ';     T=', num2str(T), ';     H=', num2str(H) ];
            disp(str);
            ISTATUS.Nsng = ISTATUS.Nsng + 1;
            ConsecutiveSng = ConsecutiveSng + 1;
            if ( ConsecutiveSng >= 6 )
                return; % Failure
            end
            H = 0.5*H;
            SkipJac = false;
            SkipLU = false;
            Reject = true;
        end
    end
    
return; 
