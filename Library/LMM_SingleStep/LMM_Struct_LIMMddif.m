function [lmms] = LMM_Struct_LIMMddif()
% Builds a LMM struct for a variable-stepsize LIMM method.

 % Precompute/preload coefficient matrices if possible.
 lmms.coefficients = @LIMMddif_Coefficients;
 
 % Compute the results for 1 step of the method.
 lmms.onestep = @LIMMddif_Onestep;

 % Initialize LIMM internal state.
 lmms.stateInit = @LIMMddif_stateInit;
 
 % Advance LIMM internal state by 1 step.
 lmms.stateAdvance = @LIMMddif_stateAdvance;
 
 % Update timestep size
 lmms.stateUpdateH = @LIMMddif_stateUpdateH;
 
end

function [state] = LIMMddif_stateInit(MaxOrder, NVAR, Y0, F0, H)
% Initialize BDF internal state and return the state object.

  state.NVAR = NVAR;
  state.K = MaxOrder + 3; % 0-K+2
  state.Ydif = [Y0, zeros(NVAR,state.K-1)];
  state.Fdif = [F0, zeros(NVAR, state.K-1)];
  state.H = [H, zeros(1,state.K-1)];
  state.YnewDif = zeros(NVAR,1);

end

function [state] = LIMMddif_stateAdvance(state, Hnew, Ynew, Fnew, Order)
% Advance LIMM internal state by one step and return the state object.

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
  
  state.Fdif(:,Order+1) = (Fnew - state.Fdif(:,1))/prod(Hsums(1:Order)) - state.Fdif(:,2:Order) * (1 ./ cumprod(Hsums(2:Order), 'reverse'))';
  for j = Order:-1:1
      state.Fdif(:,j) = state.Fdif(:,j) + Hsums(j) * state.Fdif(:,j+1);
  end
  
  % debug
  %assert(norm(state.Fdif(:,1) - Fnew)/norm(Fnew) < 100*eps, ['Assertion failed. norm(Fnew - Fdif(:,1)) = ', num2str(norm(Fnew - state.Fdif(:,1))), ', norm(Fdif(:,1)) = ', num2str(norm(state.Fdif(:,1))), '.']);

end

function [state] = LIMMddif_stateUpdateH(state, Hnew)
% Update current stepsize H.

  state.H(1) = Hnew;

end

function [c] = LIMMddif_cell_eval(w, coeffF)
% Evaluate a function with multple arguments from a single array.

    wc = num2cell(w);
    c = coeffF(wc{:});
    
end

function [coeffs] = LIMMddif_Coefficients()
% Preload all coefficient functions with direct expressions

    coeffs.maxOrder = 5;
    coeffs.c = {};
    coeffs.d = {};
    coeffs.e = {};
    coeffs.gamma = {};
    coeffs.err = {};
    for k = 1:coeffs.maxOrder
        [c, d, e, gamma, err] = limm_expr_ddif_o16(k);
        coeffs.c{k} = @(w)LIMMddif_cell_eval(w, c);
        coeffs.d{k} = @(w)LIMMddif_cell_eval(w, d);
        coeffs.e{k} = @(w)LIMMddif_cell_eval(w, e);
        coeffs.gamma{k} = @(w)LIMMddif_cell_eval(w, gamma);
        coeffs.err{k} = @(w)LIMMddif_cell_eval(w, err);
    end
    coeffs.getCoeffs = @LIMMddif_Coeffs;
    
end

function [coeffs] = LIMMddif_Coefficients_T()
% Preload all coefficient functions using Tensor form

    coeffs.maxOrder = 6;
    coeffs.c = {};
    coeffs.d = {};
    coeffs.e = {};
    coeffs.gamma = {};
    coeffs.err = {};
    for k = 1:coeffs.maxOrder
        [cstr, dstr, estr, gstr, rstr] = limm_ddif_o16(k);
        [c, d, e, gamma, err] = limm_tensor_ddif(cstr, dstr, estr, gstr, rstr);
        coeffs.c{k} = @(w)limm_tensor_eval(w, c);
        coeffs.d{k} = @(w)limm_tensor_eval(w, d);
        coeffs.e{k} = @(w)limm_tensor_eval(w, e);
        coeffs.gamma{k} = @(w)limm_tensor_eval(w, gamma);
        coeffs.err{k}  = @(w)limm_tensor_eval(w, err);
    end
    coeffs.getCoeffs = @LIMMddif_Coeffs;
    
end

function [c, d, e, gamma, err] = LIMMddif_Coeffs(Coeffs, Order, H, k)
% Returns the [c, d, e, gamma, err] coefficients for LIMM with the current steps
% for order p.

    w = H(2:end) ./ H(1);
    
    c = zeros(k,1);
    d = zeros(k,1);
    e = zeros(k,1);
    %gamma = 0;
    err = ones(1,3);
    
    c(1:Order) = Coeffs.c{Order}(w(1:Order-1))';
    d(1:Order) = Coeffs.d{Order}(w(1:Order-1))';
    e(1:Order) = Coeffs.e{Order}(w(1:Order-1))';
    gamma = Coeffs.gamma{Order}(w(1:Order-1));
    err(2) = max(abs(Coeffs.err{Order}(w(1:Order-1))));
    
    % error constants for higher and lower order methods
    if Order > 1
        err(1) = max(abs(Coeffs.err{Order-1}(w(1:Order-2))));
    end
    if Order < Coeffs.maxOrder
        err(3) = max(abs(Coeffs.err{Order+1}(w(1:Order))));
    end
    
end

function [Ynew, YE, ELO, H, RejectStep, State, ISTATUS] = LIMMddif_Onestep(Order, H, Y, F, Jacobian, dFdT, T, State, Coefficients, OdeFunction, OPTIONS, ISTATUS)

    % debug
    assert(H == State.H(1), ['H = ', num2str(H), ', State.H = ', num2str(State.H)]);
    
    %disp(['norm(H*F - Ydif(:,2))/norm(H*F) = ', num2str(norm(State.H(2)*F - State.Ydif(:,2))/norm(State.H(2)*F))])

    [c, d, e, gamma, err] = Coefficients.getCoeffs(Coefficients, Order, State.H, State.K);
    
    RejectStep = false;
    NVAR = State.NVAR;
    
    %keyboard
    Hp = H.^[0:State.K-1]';
    
    if (~OPTIONS.MatrixFree)
        RHS = -(State.Ydif * (c .* Hp)) + H * Jacobian * (State.Ydif * (d .* Hp)) + H * (State.Fdif * (e .* Hp));
        LHS = speye(NVAR) - H * gamma * Jacobian;
        Ynew = LHS\RHS;
        ISTATUS.Nsol = ISTATUS.Nsol + 1;
    else
        RHS = -(State.Ydif * (c .* Hp)) + H * Jacobian(State.Ydif * (d .* Hp)) + H * (State.Fdif * (e .* Hp));
        LHS = @(v)(v - H * gamma * Jacobian(v));
        [Ynew, ISTATUS] = MatrixFreeSolve(LHS, RHS, OPTIONS, ISTATUS);
    end
    
    % Compute error estimates
    errScale =  H.^[Order-1, Order, Order+1] .* factorial([Order-1, Order, Order+1]);
    %errScale =  H.^[Order, Order+1, Order+2] .* factorial([Order, Order+1, Order+2]);
    if ISTATUS.Nacc < 1
        State.YnewDif = (Ynew - Y)/(2*H*H) - F/(2*H);
        if ~OPTIONS.MatrixFree
            YE = [Ynew, State.YnewDif, (State.YnewDif - Jacobian*F)/(3*H)] .* err .* errScale;
        else
            YE = [Ynew, State.YnewDif, (State.YnewDif - Jacobian(F))/(3*H)] .* err .* errScale;
        end
    else
        Hsums = cumsum(State.H);
        State.YnewDif = (Ynew - Y)/prod(Hsums(1:Order+1)) - State.Ydif(:,2:Order+1) * (1 ./ cumprod(Hsums(2:Order+1), 'reverse'))';
        YE = [Hsums(Order+1) * State.YnewDif + State.Ydif(:,Order+1), State.YnewDif, (State.YnewDif - State.Ydif(:,Order+2))/Hsums(Order+2)] .* err .* errScale;
    end
    ELO = [max(1,Order-1), Order, Order+1];
    %ELO = 0:2 + Order;

end




