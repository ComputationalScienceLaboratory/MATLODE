%% ERK_ADJ1_DiscreteIntegrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ERK_ADJ1_DiscreteIntegrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )
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
% |Tout|: Time vector
%
% |Yout|: State vector
%
% |Lambda|: Sensitivity matrix
%
% |Mu|: Mu term
%
% |ISTATUS|: Integer statistics 
%
% |RSTATUS|: Real statistics
%
% |Ierr|: Error flag
%
%% Description
% Explicit Runge-Kutta discrete 1 adjoint core method.
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
function [ Tout, Yout, Lambda, Mu, ISTATUS, RSTATUS, Ierr ] = ERK_ADJ1_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag, lambda_tf, mu_tf )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Lambda = lambda_tf;
    
    % We are forcing the user to have a dim(Mu) to be NADJ x NP, but the
    % implementation expects NP x NADJ.
    Mu = transpose(mu_tf);
    
    Yout = Lambda;
    Tout = 0;
    TYindex = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    U = zeros(NVAR,OPTIONS.NADJ,Coefficient.NStage);
    V = zeros(OPTIONS.NP,OPTIONS.NADJ,Coefficient.NStage);
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( stack_ptr > 0 )
        % Recover checkpoints for stage vales and vectors
        [ T, H, Y, Z, stack_ptr ] = erkPop( NVAR, Coefficient.NStage, stack_ptr );
        Yout(:,TYindex) = Y;
        Tout(TYindex,1) = T;
        
        TYindex = TYindex + 1;
        
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
        for istage=Coefficient.NStage:-1:1
            % Jacobian of the current stage solution
            fjac = OPTIONS.Jacobian( T+rkC(istage)*H, Z(:,istage) );
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            if ( adjQuadFlag )
                WY = OPTIONS.DRDY(T+rkC(istage)*H, Z(:,istage) );
            end
            fpjac = OPTIONS.Jacp( T+rkC(istage)*H, Z(:,istage) );
            if ( adjQuadFlag )
                WP = OPTIONS.DRDP( T+rkC(istage)*H, Z(:,istage) );
                
                % This is inefficient and should be accounted for in a
                % future release
                WP = transpose(WP);
            end
            
            for iadj=1:OPTIONS.NADJ
                TMP = zeros(NVAR,1);
                
                if ( istage < Coefficient.NStage )
                    for j=istage+1:Coefficient.NStage
                        TMP = TMP + H*rkA(j,istage)*U(:,iadj,j);
                    end
                end
                TMP = TMP + H*rkB(istage)*Lambda(:,iadj);
                
                U(:,iadj,istage) = transpose(fjac)*TMP;
                V(:,iadj,istage) = transpose(fpjac)*TMP;
                
                if ( adjQuadFlag )
                    U(:,iadj,istage) = U(:,iadj,istage) + H*rkB(istage)*WY(:, iadj);
                    V(:,iadj,istage) = V(:,iadj,istage) + H*rkB(istage)*WP(:, iadj);
                end
                
            end % iadj
            
        end % istage
        
        % Update adjoint solution
        for istage=1:Coefficient.NStage
            Lambda = Lambda + U(:,:,istage);
            Mu = Mu + V(:,:,istage);
        end
        
    end % time loop
    
    % Successful return
    Ierr = 1;
    

return;

