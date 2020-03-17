function [lmms] = LMM_Struct_BDFddif()
% Builds a LMM struct for a variable-stepsize BDF method.

 % Precompute/preload coefficient matrices if possible.
 lmms.coefficients = @BDFddif_Coefficients;
 
 % Compute the results for 1 step of the method.
 lmms.onestep = @BDFddif_Onestep;

 % Initialize BDF internal state.
 lmms.stateInit = @BDFddif_stateInit;
 
 % Advance BDF internal state by 1 step.
 lmms.stateAdvance = @BDFddif_stateAdvance;
 
 % Update timestep size
 lmms.stateUpdateH = @BDFddif_stateUpdateH;
 
end

function [state] = BDFddif_stateInit(MaxOrder, NVAR, Y0, F0, H)
% Initialize BDF internal state and return the state object.

  state.NVAR = NVAR;
  state.K = MaxOrder + 3; % 0-K+2
  state.Ydif = [Y0, F0, zeros(NVAR,state.K-2)]; % Y = Ydif0, F ~ Ydif1
  state.H = H * ones(1,state.K); % w1 = 1 for first step
  state.YnewDif = zeros(NVAR,1);

end

function [state] = BDFddif_stateAdvance(state, Hnew, Ynew, ~, Order)
% Advance BDF internal state by one step and return the state object.

  Hsums = cumsum(state.H);

  % cycle in new H
  for j = state.K:-1:2
      state.H(j) = state.H(j-1);
  end
  state.H(1) = Hnew;
  
  % cycle in new Y
  state.Ydif(:,Order+3) = state.YnewDif - Hsums(Order+2) * state.Ydif(:,Order+2); % Order+2 dif approx
  state.Ydif(:,Order+2) = state.YnewDif;                         % Order+1 dif approx
  
  state.Ydif(:,Order+1) = (Ynew - state.Ydif(:,1))/prod(Hsums(1:Order)) - state.Ydif(:,2:Order) * (1 ./ cumprod(Hsums(2:Order), 'reverse'))';   % Order dif exact
  for j = Order:-1:1
      state.Ydif(:,j) = state.Ydif(:,j) + Hsums(j) * state.Ydif(:,j+1);
  end
  
  % debug
  %assert(norm(state.Ydif(:,1) - Ynew)/norm(Ynew) < 100*eps, ['Assertion failed. norm(Ynew - Ydif(:,1)) = ', num2str(norm(Ynew - state.Ydif(:,1))), ', norm(Ydif(:,1)) = ', num2str(norm(state.Ydif(:,1))), '.']);
  
end

function [state] = BDFddif_stateUpdateH(state, Hnew)
% Update current stepsize H.

  state.H(1) = Hnew;

end

function [coeffs] = BDFddif_Coefficients()
% Preload all coefficient functions with direct expressions

    coeffs.maxOrder = 5;
    coeffs.ahat = {};
    coeffs.alpha = {};
    coeffs.bhat = {};
    coeffs.beta = {};
    coeffs.err = {};
    
    for k = 1:coeffs.maxOrder
        [alpha, ahat, beta, bhat, err] = bdf_ddif_expr_o16(k);
        coeffs.alpha{k} = alpha;
        coeffs.ahat{k} = ahat;
        coeffs.beta{k} = beta;
        coeffs.bhat{k} = bhat;
        coeffs.err{k} = err;
    end
    
    coeffs.getCoeffs = @BDFddif_Coeffs;
    
end

function [alpha, ahat, beta, bhat, err] = BDFddif_Coeffs(Coeffs, Order, H, k)
% Returns the [alpha, ahat, err] coefficients for BDF with the current steps
% for order p.

    w = H(2:end) ./ H(1);
    
    err = ones(1,3);
    
    ahat = Coeffs.ahat{Order}(w(1:Order));
    alpha = Coeffs.alpha{Order}(w(1:Order))';
    bhat = Coeffs.bhat{Order}(w(1:Order));
    beta = Coeffs.beta{Order}(w(1:Order))';
    err(2) = Coeffs.err{Order}(w(1:Order-1));
    
    % error constants for higher and lower order methods
    if Order > 1
        err(1) = Coeffs.err{Order-1}(w(1:Order-2));
    end
    if Order < Coeffs.maxOrder
        err(3) = Coeffs.err{Order+1}(w(1:Order));
    end
    
    err = abs(err);
    
end

function [Ynew, YE, ELO, H, RejectStep, State, ISTATUS] = BDFddif_Onestep(Order, H, Y, F, Jacobian, dFdT, T, State, Coefficients, OdeFunction, OPTIONS, ISTATUS)

    % debug
    assert(H == State.H(1), ['H = ', num2str(H), ', State.H = ', num2str(State.H)]);
    
    %disp(['norm(H*F - Ydif(:,2))/norm(H*F) = ', num2str(norm(State.H(2)*F - State.Ydif(:,2))/norm(State.H(2)*F))])
    
    RejectStep = false;
    SkipLU = false;
    LHS = struct;
    Fac = 0.5;
    shrinkCount = 0;
    
    NVAR = State.NVAR;
    normY = norm(Y);
    
    NewtonDone = false;
    while ~NewtonDone
        
        [alpha, ahat, beta, bhat, err] = Coefficients.getCoeffs(Coefficients, Order, State.H, State.K);
        
        assert(size(alpha,1) == Order);
        assert(size(beta,1) == Order);
        
        Hp = H.^[1:Order]';
        RHSfixed = State.Ydif(:,2:Order+1) * ((Hp ./ H^(Order+1)) .* alpha);
        Ypred = Y + State.Ydif(:,2:Order+1) * (Hp .* beta);
        Ynew = Ypred;
        State.YnewDif = zeros(NVAR,1);
        
        %debug
        assert(size(RHSfixed,1) == NVAR);
        assert(size(Ypred,1) == NVAR);
        
        % compute LHS matrix LU decomposition
        if (~SkipLU)
            if (~OPTIONS.MatrixFree)
                if issparse(Jacobian)
                    J = speye(NVAR) - H * (bhat/ahat) * Jacobian;
                    [LHS.L, LHS.U, LHS.P, LHS.Q, LHS.R] = lu(J);
                else
                    J = eye(NVAR) - H * (bhat/ahat) * Jacobian;
                    [LHS.L, LHS.U, LHS.p] = lu(J, 'vector');
                end
            else % for MatrixFree operation, build J*v function handle
                LHS.jv = @(v)(v - H * (bhat/ahat) * Jacobian(v));
            end
            SkipLU = true;
        end
        
        % Begin Newton iterations
        for iter = 1:OPTIONS.NewtonMaxIt

            Fnew = OdeFunction(T+H, Ynew);
            
            RHS = 1/(H^Order * ahat) * Fnew - (State.YnewDif + RHSfixed);
            
            if (~OPTIONS.MatrixFree)
                if issparse(Jacobian)
                    dYdif = LHS.Q * (LHS.U \ (LHS.L \ (LHS.P * (LHS.R \ RHS))));
                else
                    dYdif = LHS.U \ (LHS.L \ RHS(LHS.p));
                end
                ISTATUS.Nsol = ISTATUS.Nsol + 1;
            else
                [dYdif, ISTATUS] = MatrixFreeSolve(LHS.jv, RHS, OPTIONS, ISTATUS);
            end
            
            State.YnewDif = State.YnewDif + dYdif;
            Ynew = Ypred + H^(Order+1) * bhat * State.YnewDif;
            
            newnorm = H^(Order+1) * bhat * norm(dYdif);
            
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
            
        end % NewtonIter
        
        if ( ~NewtonDone )
            if (shrinkCount > 5)
                warning('Newton failed to converge after more than 5 attempts.');
                RejectStep = true;
                break;
            else
                %disp(['Shrinking stepsize... Fac = ', num2str(Fac), '.']);
                shrinkCount = shrinkCount + 1;
                State.H(1) = max(OPTIONS.Hmin, min(Fac*H(1), OPTIONS.Hmax));
                H = State.H(1);
            
                % Could also have logic to recompute Jacobian here...
            end
        end
        
    end % NewtonDone
    
    % Compute error estimates
    %errScale =  err .* H.^[Order-1, Order, Order+1] .* factorial([Order-1, Order, Order+1]);
    errScale = err .* H.^[Order, Order+1, Order+2] .* factorial([Order, Order+1, Order+2]); % This one SHOULD be right, but....
    if ISTATUS.Nacc < 1 % No previous steps, so assume Order == 1
        if ~OPTIONS.MatrixFree
            YE = [ones(size(Ynew)), errScale(2) * State.YnewDif, errScale(3) * (State.YnewDif - Jacobian*F)/(3*H)];
        else
            YE = [ones(size(Ynew)), errScale(2) * State.YnewDif, errScale(3) * (State.YnewDif - Jacobian(F))/(3*H)];
        end
    else
        Hsums = cumsum(State.H);
        YE = [errScale(1) * (Hsums(Order+1) * State.YnewDif + State.Ydif(:,Order+1)), errScale(2) * State.YnewDif, errScale(3) * (State.YnewDif - State.Ydif(:,Order+2))/Hsums(Order+2)];
    end
    ELO = [max(1,Order-1), Order, Order+1];
    %ELO = 0:2 + Order;

end




