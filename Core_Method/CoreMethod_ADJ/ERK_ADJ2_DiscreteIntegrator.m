%% ERK_ADJ2_DiscreteIntegrator
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
function [ Lambda, ISTATUS, RSTATUS, Ierr ] = ERK_ADJ2_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC

    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Required User Input
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( isempty(OPTIONS.Lambda) ) 
        error('User defined option parameter Lambda is required.');
    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    U = zeros(NVAR,OPTIONS.NADJ,Coefficient.NStage);
    
    Lambda = OPTIONS.Lambda;
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( stack_ptr > 0 ) % Tloop
        % Recover checkpoints for stage vales and vectors
        [ T, H, ~, Z, stack_ptr ] = erkPop( NVAR, Coefficient.NStage, stack_ptr, OPTIONS );
                
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
        for istage = Coefficient.NStage:-1:1
            % Jacobian of the current stage solution
            fjac = OPTIONS.Jacobian(T+rkC(istage)*H,Z(:,istage));
            ISTATUS.Njac = ISTATUS.Njac + 1;

            if ( adjQuadFlag ) 
                WY = OPTIONS.DRDY(T+rkC(istage)*H,Z(:,istage));
            end
                        
            for iadj=1:OPTIONS.NADJ % adj
                TMP = zeros(NVAR,1);
                if ( istage < Coefficient.NStage )
                    for j=istage+1:Coefficient.NStage
                        TMP = TMP + H*rkA(j,istage)*U(:,iadj,j); % TMP < h*A*U
                    end
                end

                TMP = TMP + H*rkB(istage)*Lambda(:,iadj);
                
                U(:,iadj,istage) = transpose(fjac)*TMP;

                if ( adjQuadFlag )
                    U(:,iadj,istage) = U(:,iadj,istage) + H*rkB(istage)*WY(:,iadj);
                end
            end % adj
        end % stages

        % Update adjoint solution
        % Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
        for istage=1:Coefficient.NStage
            for iadj=1:OPTIONS.NADJ
                Lambda(:,iadj) = Lambda(:,iadj) + U(:,iadj,istage);  % <----- original line
            end
        end
    end % Tloop
    
    % Successful return
    Ierr = 1;

return;

function [ T, H, Y, Z, stack_ptr ] = erkPop( NVAR, rkS, stack_ptr, OPTIONS )
    global chk_H chk_T chk_Y chk_Z

    if ( stack_ptr <= 0 )
        error( 'Pop failed: empty buffer' );
    end
    H = chk_H( stack_ptr );
    T = chk_T( stack_ptr );
    
     % for debugging
     if ( OPTIONS.displaySteps == true )
        str = ['Backtracking. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
        disp(str);
     end    
    
    Y(1:NVAR) = chk_Y( 1:NVAR, stack_ptr );
    Z(1:NVAR, 1:rkS) = chk_Z(1:NVAR, 1:rkS, stack_ptr );
    %Z(1:rkS,1:NVAR) = chk_Z(:, :, stack_ptr );
    
    stack_ptr = stack_ptr - 1;

return;
