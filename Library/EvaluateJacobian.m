function [ fjac, ISTATUS ] = EvaluateJacobian(T, Y, F, OdeFunction, OPTIONS, ISTATUS)
%EvaluateJacobian Computes the Jacobian function J(T, Y)
% Or, if MatrixFree, builds J(T,Y)*v Jacobian-vector product handle
    if ( ~OPTIONS.MatrixFree )
        fjac = OPTIONS.Jacobian(T, Y);
        ISTATUS.Njac = ISTATUS.Njac + 1;
    else
        if ( ~isempty( OPTIONS.Jacobian ) )
            if ( nargin( OPTIONS.Jacobian ) == 3 ) % Jacobian-vector product function supplied
                fjac = @(vee)OPTIONS.Jacobian(T, Y, vee);
            elseif ( nargin( OPTIONS.Jacobian ) == 2 )
                Jac = OPTIONS.Jacobian(T, Y);
                %ISTATUS.Njac = ISTATUS.Njac + 1;
                fjac = @(vee)(Jac*vee);
            else
                error('The Jacobian function should require only 2 or 3 arguments.')
            end
        else
            if ( isempty( F ) )
                F = OdeFunction(T, Y);
                ISTATUS.Nfun = ISTATUS.Nfun+1;
            end
            normy = norm(Y);
            fjac = @(v)Mat_Free_Jac(T, Y, v, OdeFunction, F, normy);
        end
    end

end

