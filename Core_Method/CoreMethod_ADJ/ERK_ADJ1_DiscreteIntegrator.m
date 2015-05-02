%% ERK_ADJ1_DiscreteIntegrator
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
function [ Tout, Yout, Lambda, Mu, ISTATUS, RSTATUS, Ierr ] = ERK_ADJ1_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Lambda = OPTIONS.Lambda;
    Yout = Lambda;
    Tout = 0;
    TYindex = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    U = zeros(NVAR,OPTIONS.NADJ,Coefficient.NStage);
    V = zeros(OPTIONS.NP,OPTIONS.NADJ,Coefficient.NStage);
    
    Mu = OPTIONS.Mu();

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
                    U(:,iadj,istage) = U(:,iadj,istage) + H*rkB(istage)*WY(:,iadj);
                    V(:,iadj,istage) = V(:,iadj,istage) + H*rkB(istage)*WP(:,iadj);
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

