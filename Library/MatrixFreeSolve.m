function [X, ISTATUS, gmresFlag, singFlag] = MatrixFreeSolve(LHS, RHS, OPTIONS, ISTATUS)
% wrap gmres and it's warnings
    [X, gmresFlag, ~, iter] = gmres(LHS, RHS, OPTIONS.GMRES_Restart, OPTIONS.GMRES_TOL, min(OPTIONS.GMRES_MaxIt,length(RHS)), OPTIONS.GMRES_P);
    ISTATUS.Nsol = ISTATUS.Nsol + 1;

    if ( ~isempty(OPTIONS.GMRES_Restart) )
        vecCount = iter(2) + (OPTIONS.GMRES_Restart - 1)*iter(1);
    else
        vecCount = iter(2);
    end
    ISTATUS.Njac =  ISTATUS.Njac + vecCount;

    singFlag = 0;
    if( gmresFlag ~= 0 )
        resvec = abs(LHS(X) - RHS);
        scalar = OPTIONS.AbsTol + OPTIONS.RelTol.*abs(RHS);
        if (norm(resvec./scalar) > sqrt(length(RHS)))
            singFlag = 1;
            ISTATUS.Nsng = ISTATUS.Nsng + 1;
            switch(gmresFlag)
                case 1
                    warning('GMRES: iterated MAXIT times but did not converge');
                case 2
                    warning('GMRES: preconditioner M was ill-conditioned');
                case 3
                    warning('GMRES: stagnated (two consecutive iterates were the same)');
            end
        else
            gmresFlag = 0;
        end
    end
end
