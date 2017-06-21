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
function [ Tout, Yout, Lambda, Mu, ISTATUS, RSTATUS, Ierr ] = ERK_ADJ_DiscreteIntegrator(Lambda_tf,Mu_tf, OPTIONS, Coefficient, stack_ptr, adjQuadFlag, adjMuFlag )

% adjMuFlag says if Mu is calculated, if so, Jacp is necessary and DRDP is
% required if a quadrature variable z=integral r is provided

% adjquadFlag says if Lambda  is calculated with DRDY, where a the sensitivity of the quadrature
% is provided



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global rkA rkB rkC

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tout = 0;

TYindex = 1;
Lambda=transpose(Lambda_tf);
Mu=transpose(Mu_tf);
%     Yout = Lambda;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

U = zeros(OPTIONS.NVAR,OPTIONS.NADJ,Coefficient.NStage);
if adjMuFlag
    V = zeros(OPTIONS.NP,OPTIONS.NADJ,Coefficient.NStage);
end


ISTATUS = ISTATUS_Struct('default');
RSTATUS = RSTATUS_Struct('default');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while ( stack_ptr > 0 )
    % Recover checkpoints for stage vales and vectors
    [ T, H, Y, Z, stack_ptr ] = erkPop( OPTIONS.NVAR, Coefficient.NStage, stack_ptr );
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = T;
    
    TYindex = TYindex + 1;
    
    ISTATUS.Nstp = ISTATUS.Nstp + 1;
    
    for istage=Coefficient.NStage:-1:1
        % Jacobian of the current stage solution
        fjac = OPTIONS.Jacobian( T+rkC(istage)*H, Z(:,istage) );
        ISTATUS.Njac = ISTATUS.Njac + 1;
        
        if  adjQuadFlag
            WY = OPTIONS.DRDY(T+rkC(istage)*H, Z(:,istage) );
            WY=transpose(WY);
        end
        
        if  adjMuFlag %  Jacp and DRDP are used in the Mu calculation
            fpjac = OPTIONS.Jacp( T+rkC(istage)*H, Z(:,istage) );
            if  adjQuadFlag
            WP = OPTIONS.DRDP( T+rkC(istage)*H, Z(:,istage) );
            WP=transpose(WP);
            end
        end
        
        for iadj=1:OPTIONS.NADJ
            
            TMP = zeros(OPTIONS.NVAR,1);            
            if ( istage < Coefficient.NStage )
                for j=istage+1:Coefficient.NStage
                    TMP = TMP + H*rkA(j,istage)*U(:,iadj,j);
                end
            end
            TMP = TMP + H*rkB(istage)*Lambda(:,iadj);
            U(:,iadj,istage) = transpose(fjac)*TMP;
            
            if  adjQuadFlag
                U(:,iadj,istage) = U(:,iadj,istage) + H*rkB(istage)*WY(:,iadj);
            end
            if adjMuFlag
                V(:,iadj,istage) = transpose(fpjac)*TMP;
                if adjQuadFlag
                V(:,iadj,istage) = V(:,iadj,istage) + H*rkB(istage)*WP(:,iadj);
                end
            end
        end % iadj
        
    end % istage
    
    % Update adjoint solution
    for istage=1:Coefficient.NStage
        Lambda = Lambda + U(:,:,istage);
        if adjMuFlag
        Mu = Mu + V(:,:,istage);
        end
    end
    
end % time loop

% Successful return
Ierr = 1;
Mu=transpose(Mu);
Lambda=transpose(Lambda);

return;

