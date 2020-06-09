function [lmms] = LMM_Struct_LIMM()
% Builds a LMM struct for a variable-stepsize LIMM method.

 % Precompute/preload coefficient matrices if possible.
 lmms.coefficients = @LIMM_Coefficients;
 
 % Compute the results for 1 step of the method.
 lmms.onestep = @LIMM_Onestep;

 % Initialize LIMM internal state.
 lmms.stateInit = @LIMM_stateInit;
 
 % Advance LIMM internal state by 1 step.
 lmms.stateAdvance = @LIMM_stateAdvance;
 
 % Update timestep size
 lmms.stateUpdateH = @LIMM_stateUpdateH;
 
end

function [state] = LIMM_stateInit(MaxOrder, NVAR, Y0, F0, H)
% Initialize LIMM internal state and return the state object.

  state.K = MaxOrder + 1;
  state.Y = [zeros(NVAR, state.K-1), Y0];
  state.F = [zeros(NVAR, state.K-1), F0];
  state.H = [zeros(1,state.K-1), H];
  state.IDX = 1:state.K; % Y(:,IDX(k)) is always the current Y, and accepted Ynew are stored in Y(IDX(1)).

end

function [state] = LIMM_stateAdvance(state, Hnew, Ynew, Fnew, ~)
% Advance LIMM internal state by one step and return the state object.

  state.Y(:,state.IDX(1)) = Ynew;
  state.F(:,state.IDX(1)) = Fnew;
  state.H(state.IDX(1)) = Hnew;
  
  tmp = state.IDX(1);
  for j = 1:state.K-1
      state.IDX(j) = state.IDX(j+1);
  end
  state.IDX(state.K) = tmp;

end

function [state] = LIMM_stateUpdateH(state, Hnew)
% Update stepsize H after rejected step.

  state.H(state.IDX(state.K)) = Hnew;

end

function [c] = LIMM_cell_eval(w, coeffF)
% Evaluate a function with multple arguments from a single array.

    wc = num2cell(w);
    c = coeffF(wc{:});
    
end

function [coeffs] = LIMM_Coefficients()
% Preload all coefficient functions with direct expressions

    coeffs.maxOrder = 6;
    coeffs.alpha = {};
    coeffs.beta  = {};
    coeffs.gamma = {};
    for k = 1:coeffs.maxOrder
        [alpha, beta, gamma] = limm_expr_o16(k);
        coeffs.alpha{k} = @(w)LIMM_cell_eval(w, alpha);
        coeffs.beta{k}  = @(w)LIMM_cell_eval(w, beta);
        coeffs.gamma{k} = @(w)LIMM_cell_eval(w, gamma);
    end
    coeffs.getCoeffs = @LIMM_Coeffs;
    
end

function [coeffs] = LIMM_Coefficients_T()
% Preload all coefficient functions using Tensor form

    coeffs.maxOrder = 6;
    coeffs.alpha = {};
    coeffs.beta  = {};
    coeffs.gamma = {};
    for k = 1:coeffs.maxOrder
        [astr, bstr, gstr] = limm_o16(k);
        [alpha, beta, gamma] = limm_tensor(astr, bstr, gstr);
        coeffs.alpha{k} = @(w)limm_tensor_eval(w, alpha);
        coeffs.beta{k}  = @(w)limm_tensor_eval(w, beta);
        coeffs.gamma{k} = @(w)limm_tensor_eval(w, gamma);
    end
    coeffs.getCoeffs = @LIMM_Coeffs;
    
end

function [alphas, betas, gammas, gam] = LIMM_Coeffs(Coeffs, Order, H, k)
% Returns the [alpha, beta, gamma] coefficients for LIMM with the current steps
% for orders: O-1,O,O+1. Always uses the same coefficient family (o12,o26).

    w = H(1:end-1) ./ H(end);
    
    if (Order == 1)
        apm1 = []; bpm1 = []; gpm1 = [];
    else
        apm1 = Coeffs.alpha{Order-1}(w(end-Order+3:end));
        bpm1 = Coeffs.beta{Order-1}(w(end-Order+3:end));
        gpm1 = Coeffs.gamma{Order-1}(w(end-Order+3:end));
    end
    ap = Coeffs.alpha{Order}(w(end-Order+2:end));
    bp = Coeffs.beta{Order}(w(end-Order+2:end));
    gp = Coeffs.gamma{Order}(w(end-Order+2:end));
    if (Order < Coeffs.maxOrder)
        app1 = Coeffs.alpha{Order+1}(w(end-Order+1:end));
        bpp1 = Coeffs.beta{Order+1}(w(end-Order+1:end));
        gpp1 = Coeffs.gamma{Order+1}(w(end-Order+1:end));
    else
        app1 = []; bpp1 = []; gpp1 = [];
    end
    
    %keyboard
    
    alphas = zeros(k, 3);
    alphas(end-length(apm1)+1:end,1) = apm1';
    alphas(end-length(ap)+1:end,2)   = ap';
    alphas(end-length(app1)+1:end,3) = app1';
    
    betas = zeros(k, 3);
    betas(end-length(bpm1)+1:end,1) = bpm1';
    betas(end-length(bp)+1:end,2)   = bp';
    betas(end-length(bpp1)+1:end,3) = bpp1';
    
    gam = zeros(1,3);
    if (Order > 1)
        gam(1) = gpm1(end);
    end
    gam(2) = gp(end);
    if (Order < Coeffs.maxOrder)
        gam(3) = gpp1(end);
    end
    gammas = zeros(k, 3);
    gammas(end-length(gpm1(1:end-1))+1:end,1) = gpm1(1:end-1)';
    gammas(end-length(gp(1:end-1))+1:end,2)   = gp(1:end-1)';
    gammas(end-length(gpp1(1:end-1))+1:end,3) = gpp1(1:end-1)';
    
end

function [Yp, YE, ELO, H, RejectStep, State, ISTATUS] = LIMM_Onestep(Order, H, Y, F, Jacobian, dFdT, T, State, Coefficients, OdeFunction, OPTIONS, ISTATUS)

% Debug
assert(H == State.H(State.IDX(State.K)),['H = ', num2str(H), ' State.H = ', num2str(State.H(State.IDX)), ' diff H = ', num2str(H - State.H(State.IDX(State.K)))])
assert(sum(Y - State.Y(:,State.IDX(State.K))) < 20*eps)

    % shorthand
    IDX = State.IDX;
    k = State.K;

    [alpha, beta, gamma, gam] = Coefficients.getCoeffs(Coefficients, Order, State.H(IDX), k);
    
    methods = 1:3;
    if (Order == 1)
        methods = methods(2:end);
    end
    if (Order + 1 > nnz(State.H))
        methods = methods(1:end-1);
    end
    
    RejectStep = false;
    NVAR = length(Y);
    Ynew = zeros(NVAR,3);
    
    %keyboard
    
    for i = methods
        
        RHS = -(State.Y(:,IDX) * alpha(:,i)) + H * (State.F(:,IDX) * beta(:,i)) + H * Jacobian * (State.Y(:,IDX) * gamma(:,i));

        if (~OPTIONS.MatrixFree)
            LHS = speye(NVAR) - H * gam(i) * Jacobian;
            Ynew(:,i) = LHS\RHS;
            ISTATUS.Nsol = ISTATUS.Nsol + 1;
        else
            LHS = @(v)(v - H * gam(i) * Jacobian(v));
            [dY, ISTATUS] = MatrixFreeSolve(LHS, RHS, OPTIONS, ISTATUS);
            Ynew(:,i) = dY;
        end
    end
    
    Yp = Ynew(:,2);
    
    % Compute error estimates
    YE = zeros(NVAR,3);
    ELO = ones(1,3);
    
    if ( nnz(State.H) < 2 ) % for the first step only
        assert(Order == 1, 'Order is not 1 for the first step.');
        [Ypp1, failureFlag, ISTATUS] = Rosenbrock2_onestep(H, Y, F, Jacobian, T, OdeFunction, OPTIONS, ISTATUS);
        Ynew(:,3) = Ypp1;
        if (failureFlag)
            warning('Rosenbrock2 failed to converge.');
            RejectStep = true;
        end
    end
    
%     YE(:,1) = Ynew(:,3) - Ynew(:,1); ELO(1) = Order;
    YE(:,1) = Yp - Ynew(:,1); ELO(1) = Order;
    YE(:,2) = Ynew(:,3) - Yp; ELO(2) = Order+1;
    YE(:,3) = Ynew(:,3) - Yp; ELO(3) = Order+1;
%     if Order ~= 1
%         YE(:,3) = Ynew(:,3) - Ynew(:,1);
%         ELO(3) = Order;
%     else
%         YE(:,3) = Ynew(:,3) - Yp;
%         ELO(3) = Order+1;
%     end
    
%     opts = odeset('AbsTol', 100*eps, 'RelTol', 100*eps);
%     [~,Y45] = ode45(OdeFunction, [T T+H(IDX(k))], Y(:,IDX(k)), opts);
%     
%     YE = Ynew - repmat(Y45(end,:)', [1,3]);
%     ELO = [Order, Order+1, min(Order+2,6)];
    
end

function [Y2, failed, ISTATUS] = Rosenbrock2_onestep(h, y, f, fjac, T, OdeFunction, OPTIONS, ISTATUS)
% Evaluates a single step of a 2nd order Rosenbrock-W method.

    gamma = 0.282893218813449984772;
    gamma21 = -0.585786437626899969544;
%     gamma = 2.928932188134e-1;
%     gamma21 = -5.857864376269e-1;
    alpha21 = 1;
    b = 1/2;

    if ( ~OPTIONS.MatrixFree )
        lhs = speye(length(y)) - h * gamma * fjac;
        k1 = lhs\(h * f);
        f2 = OdeFunction(T + alpha21 * h, y + alpha21 * k1);
        k2 = lhs\(h * f2 + h * gamma21 * fjac * k1);
        Y2 = y + b * (k1 + k2);
        ISTATUS.Nsol = ISTATUS.Nsol + 2;
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
    else
        lhs = @(v)(v - h * gamma * fjac(v));
        rhs = h * f;
        [k1, ISTATUS] = MatrixFreeSolve(lhs, rhs, OPTIONS, ISTATUS);
        f2 = OdeFunction(T + alpha21 * h, y + alpha21 * k1);
        rhs = h * f2 + h * gamma21 * fjac(k1);
        [k2, ISTATUS] = MatrixFreeSolve(lhs, rhs, OPTIONS, ISTATUS);
        Y2 = y + b * (k1 + k2);
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        ISTATUS.Njac = ISTATUS.Njac + 1;
    end
    
    failed = false;

end


