function [lmms] = LMM_Struct_BDF()
% Builds a LMM struct for a variable-stepsize BDF method.

 % Precompute/preload coefficient matrices if possible.
 lmms.coefficients = @BDF_Coefficients;
 
 % Compute the results for 1 step of the method.
 lmms.onestep = @BDF_Onestep;
 
 % Initialize BDF internal state.
 lmms.stateInit = @BDF_stateInit;
 
 % Advance BDF internal state by 1 step.
 lmms.stateAdvance = @BDF_stateAdvance;
 
% Update timestep size
 lmms.stateUpdateH = @BDF_stateUpdateH;

 end

function [state] = BDF_stateInit(MaxOrder, NVAR, Y0, H)
% Initialize BDF internal state and return the state object.

  state.K = MaxOrder + 1;
  state.Y = [zeros(NVAR, state.K-1), Y0];
  state.H = [zeros(1,state.K-1), H];
  state.IDX = 1:state.K; % Y(IDX(k)) is always the current Y, and accepted Ynew are stored in Y(IDX(1)).

end

function [state] = BDF_stateAdvance(state, Hnew, Ynew)
% Advance BDF internal state by one step and return the state object.

  state.Y(:,state.IDX(1)) = Ynew;
  state.H(state.IDX(1)) = Hnew;
  
  tmp = state.IDX(1);
  for j = 1:state.K-1
      state.IDX(j) = state.IDX(j+1);
  end
  state.IDX(state.K) = tmp;

end

function [state] = BDF_stateUpdateH(state, Hnew)
% Update stepsize H after rejected step.

  state.H(state.IDX(state.K)) = Hnew;

end

function [coeffs] = BDF_Coefficients()
% Preload all coefficient functions

  coeffs.maxOrder = 6;
  coeffs.alpha_o = {};
  for k = 1:coeffs.maxOrder
      coeffs.alpha_o{k} = bdf_expr_o16(k);
  end
  coeffs.alpha = @BDF_Coeff_Alpha;

end

function [alphas] = BDF_Coeff_Alpha(Coeffs, Order, H, k)
% Returns the alpha coefficients for BDF of with the current steps
% for orders: O-1,O,O+1

    w = H ./ H(end);
    wc = num2cell(w);

    if (Order > 1)
        apm1 = Coeffs.alpha_o{Order-1}(wc{end-Order+2:end-1});
    else
        apm1 = [];
    end
    if (Order < Coeffs.maxOrder)
        app1 = Coeffs.alpha_o{Order+1}(wc{end-Order:end-1});
    else
        app1 = [];
    end
    ap = Coeffs.alpha_o{Order}(wc{end-Order+1:end-1});

    alphas = zeros(k+1, 3);
    alphas(end-length(apm1)+1:end,1) = apm1';
    alphas(end-length(ap)+1:end,2)   = ap';
    alphas(end-length(app1)+1:end,3) = app1';

end

function [Yp, YE, ELO, H, RejectStep, State, ISTATUS] = BDF_Onestep(Order, H, Y, F, Jacobian, dFdT, T, State, Coefficients, OdeFunction, OPTIONS, ISTATUS)

% a_(n+1) y_(n+1) + sum_i a_(n-i+1) y_(n-i+1) - h_n f(y_(n+1)) = 0

% Debug
assert(H == State.H(State.IDX(State.K)),['H = ', num2str(H), ' State.H = ', num2str(State.H(State.IDX)), ' diff H = ', num2str(H - State.H(State.IDX(State.K)))])
assert(sum(Y - State.Y(:,State.IDX(State.K))) < 20*eps)

    methods = 1:3;
    if (Order == 1)
        methods = methods(2:end);
    end
    if (Order + 1 > nnz(State.H))
        methods = methods(1:end-1);
    end
    
    % shorthand
    k = State.K;
    IDX = State.IDX;
    
    RejectStep = false;
    ShrinkStep = false;
    shrinkCount = 0;
    Fac = 0.3;
    
    normY = norm(Y);
    
    NVAR = length(Y);
    
    %keyboard
    
    NewtonDone = false;
    while ~NewtonDone

        % get coefficients for current stepsizes
        alpha = Coefficients.alpha(Coefficients, Order, State.H(IDX), k);

        Ynew = repmat(Y, [1,3]);
        
        % solve for 1 method at a time
        for j = methods

            Fnew = F;
            Jac = Jacobian;
            SkipLU = false;
            SkipJac = true;
            LHS = struct;
            
            % build Newton iteration invariant RHS
            RHSfixed = State.Y(:,IDX) * alpha(1:end-1,j);
            
            Ynew(:,j) = (1/alpha(end,j)) * (H*Fnew - RHSfixed);

            % Begin Newton iterations
            for iter = 1:OPTIONS.NewtonMaxIt

                % evaluate Jacobian
                if (~SkipJac)
                    [Jac, ISTATUS] = EvaluateJacobian(T+H, Ynew(:,j), Fnew, OPTIONS, ISTATUS);
                    SkipJac = true;
                end
                
                % compute LHS matrix LU decomposition
                if (~SkipLU)
                    if (~OPTIONS.MatrixFree)
                        if issparse(Jacobian)
                            J = alpha(end,j) * speye(NVAR) - H * Jac;
                            [LHS.L, LHS.U, LHS.P, LHS.Q, LHS.R] = lu(J);
                        else
                            J = alpha(end,j) * eye(NVAR) - H * Jac;
                            [LHS.L, LHS.U, LHS.p] = lu(J, 'vector');
                        end
                    else % for MatrixFree operation, build J*v function handle
                        LHS.jv = @(v)(alpha(end,j) * v - H * Jac(v));
                    end
                    SkipLU = true;
                end
                
                % evaluate F
                Fnew = OdeFunction(T + H, Ynew(:,j));
                
                % build RHS
                RHS = -(Ynew(:,j) * alpha(end,j) - H * Fnew + RHSfixed);

                % solve linear system
                if (~OPTIONS.MatrixFree)
                    if issparse(Jacobian)
                        dY = LHS.Q * (LHS.U \ (LHS.L \ (LHS.P * (LHS.R \ RHS))));
                    else
                        dY = LHS.U \ (LHS.L \ RHS(LHS.p));
                    end
                    ISTATUS.Nsol = ISTATUS.Nsol + 1;
                else
                    [dY, ISTATUS] = MatrixFreeSolve(LHS.jv, RHS, OPTIONS, ISTATUS);
                end

%                  % Check convergence of Newton iterations
%                 NewtonIncrement = errorNorm( NVAR, dY, SCAL );
%                 if ( iter == 1 )
%                     Theta = abs(OPTIONS.ThetaMin);
%                     NewtonRate = 2.0;
%                 else
%                     Theta = NewtonIncrement/NewtonIncrementOld;
%                     if ( Theta < 0.99 )
%                         NewtonRate = Theta/(1.0-Theta);
%                         % Predict error at the end of Newton process
%                         NewtonPredictedErr = NewtonIncrement*Theta^( OPTIONS.NewtonMaxIt - iter )/( 1.0 - Theta );
%                         if ( NewtonPredictedErr >= OPTIONS.NewtonTol )
%                             % Non-convergence of Newton: predicted error too large
%                             Qnewton = min( 10.0, NewtonPredictedErr/OPTIONS.NewtonTol );
%                             Fac = 0.8*Qnewton^(-1.0/( 1 + OPTIONS.NewtonMaxIt - iter ) );
%                             break; % NewtonLoop (confirm this)
%                         end
%                     else % Non-convergence of Newton: Theta too large
%                         break; % NewtonLoop
%                     end
%                 end
%                 NewtonIncrementOld = NewtonIncrement;

                newnorm = norm(dY);

                Ynew(:,j) = Ynew(:,j) + dY;
%                 Fnew = OdeFunction(T + H(IDX(k)), Ynew(:,j));
                
                if ( newnorm <= normY*OPTIONS.NewtonTol )
                    NewtonDone = true;
                    break;
                elseif ( iter == 1 )
                    % nothing
                elseif ( newnorm > 0.9*oldnorm )
                    % too slow
                    break;
                end
                oldnorm = newnorm;

%                 NewtonDone = ( NewtonRate*NewtonIncrement <= OPTIONS.NewtonTol );
%                 if ( NewtonDone )
%                     break; % NewtonLoop
%                 end
                
%                 SkipJac = false;
%                 SkipLU = false;
                
            end % End NewtonLoop

            if ( ~NewtonDone )
                ShrinkStep = true;
                break;
            end
        end % End method loop

        if ( shrinkCount > 5 )
            warning('Newton failed to converge after more than 5 attempts.');
            RejectStep = true;
            break;
        end
        if ( ShrinkStep )
            disp(['Im shrinking... Fac = ', num2str(Fac), '.']);
            shrinkCount = shrinkCount + 1;
            State.H(IDX(k)) = max(OPTIONS.Hmin, min(Fac*H(IDX(k)), OPTIONS.Hmax));
            H = State.H(IDX(k));
        end

    end % NewtonDone loop
    
    Yp   = Ynew(:,2);
    
    % Compute error estimates
    YE = zeros(NVAR,3);
    ELO = ones(1,3);
    if ( nnz(State.H) < 2 ) % for the first step only
        assert(Order == 1, 'Order is not 1 for the first step.');
        [Ypp1, failureFlag, ISTATUS] = Midpoint_onestep(H, Y, F, Jacobian, T, OdeFunction, OPTIONS, ISTATUS);

        YE(:,1) = Yp - Ynew(:,1);
        YE(:,2) = Ypp1 - Yp;
        YE(:,3) = Yp - Ynew(:,3);
        
        ELO(1) = 1;
        ELO(2) = 2;
        ELO(3) = 1;
        
        if (failureFlag)
            warning('Midpoint rule Newton iteration failed to converge.');
            RejectStep = true;
        end
    else
%        YE(:,1) = Ynew(:,3) - Ynew(:,1); ELO(1) = Order;
        YE(:,1) = Yp - Ynew(:,1); ELO(1) = Order;
        YE(:,2) = Ynew(:,3) - Yp; ELO(2) = Order+1;
        YE(:,3) = Ynew(:,3) - Yp; ELO(3) = Order+1;
%         if Order ~= 1
%             YE(:,3) = Ynew(:,3) - Ynew(:,1);
%             ELO(3) = Order;
%         else
%             YE(:,3) = Ynew(:,3) - Yp;
%             ELO(3) = Order+1;
%         end
    end

end

function [Y2, failed, ISTATUS] = Midpoint_onestep(h, y, f, fjac, T, OdeFunction, OPTIONS, ISTATUS)
% Evaluates a single step of the 2nd order midpoint rule.

    normY = norm(y);
    nonconv = false;
    NVAR = length(y);

    k1 = f;
    
    LHS = struct;
    if (~OPTIONS.MatrixFree)
        if issparse(fjac)
            J = speye(NVAR) - 0.5 * h * fjac;
            [LHS.L, LHS.U, LHS.P, LHS.Q, LHS.R] = lu(J);
        else
            J = eye(NVAR) - 0.5 * h * fjac;
            [LHS.L, LHS.U, LHS.p] = lu(J, 'vector');
        end
    else % for MatrixFree operation, build J*v function handle
        LHS.jv = @(v)(v - 0.5 * h * fjac(v));
    end
    
    for iter = 1:OPTIONS.NewtonMaxIt
        
        fnew = OdeFunction(T + 0.5 * h, y + 0.5 * h * k1);
        RHS = -(k1 - fnew);
        
        if ( ~OPTIONS.MatrixFree )
            if issparse(fjac)
                dK = LHS.Q * (LHS.U \ (LHS.L \ (LHS.P * (LHS.R \ RHS))));
            else
                dK = LHS.U \ (LHS.L \ RHS(LHS.p));
            end
            ISTATUS.Nsol = ISTATUS.Nsol + 1;
        else
            [dK, ISTATUS] = MatrixFreeSolve(LHS.jv, RHS, OPTIONS, ISTATUS);
        end
        
        newnorm = norm(dK);
        k1 = k1 + dK;
        
        if ( newnorm <= normY*OPTIONS.NewtonTol )
            break;
        elseif ( iter == 1 )
            % nothing
        elseif ( newnorm > 0.9*oldnorm )
            nonconv = true;
            break;
        end
        oldnorm = newnorm;
    end

    failed = nonconv;
    Y2 = y + h * k1;
    
end

