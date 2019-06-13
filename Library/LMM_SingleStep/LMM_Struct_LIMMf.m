function [lmms] = LMM_Struct_LIMMf()
% Builds a LMM struct for a variable-stepsize BDF method.

 % Precompute/preload coefficient matrices if possible.
 lmms.coefficients = @LIMMf_Coefficients;
 
 % Compute the results for 1 step of the method.
 lmms.onestep = @LIMMf_Onestep;

end


function [coeffs] = LIMMf_Coefficients()
% Preload all coefficient functions

    coeffs.maxOrder = 6;
    coeffs.alpha_o26 = {};
    coeffs.beta_o26  = {};
    coeffs.gamma_o26 = {};
    for k = 1:coeffs.maxOrder
        [astr, bstr, gstr] = limm_o26(k);
        [alpha, beta, gamma] = limm_tensor_o26(astr, bstr, gstr);
        coeffs.alpha_o26{k} = @(w)limm_tensor_eval(w, alpha);
        coeffs.beta_o26{k}  = @(w)limm_tensor_eval(w, beta);
        coeffs.gamma_o26{k} = @(w)limm_tensor_eval(w, gamma);
    end
    coeffs.transitionOrder = 2;
    coeffs.alpha_o12 = {};
    coeffs.beta_o12  = {};
    coeffs.gamma_o12 = {};
    for k = 1:coeffs.transitionOrder
        [alpha, beta, gamma] = limm_expr_o12(k);
        coeffs.alpha_o12{k} = alpha;
        coeffs.beta_o12{k}  = beta;
        coeffs.gamma_o12{k} = gamma;
    end
    coeffs.getCoeffs = @LIMMf_Coeffs;
    
end

function [alphas, betas, gammas, gam] = LIMMf_Coeffs(Coeffs, Order, H, k)
% Returns the [alpha, beta, gamma] coefficients for LIMM with the current steps
% for orders: O-1,O,O+1. Always uses the same coefficient family (o12,o26).

    w = H ./ H(end); % flip or no?
    wc = num2cell(w);
    
    if (true && Order == 1) % use family o12
        apm1 = []; bpm1 = []; gpm1 = [];
        ap   = Coeffs.alpha_o12{1}(wc{end-Order+1:end-1});
        bp   = Coeffs.beta_o12{1}(wc{end-Order+1:end-1});
        gp   = Coeffs.gamma_o12{1}(wc{end-Order+1:end-1});
        app1 = Coeffs.alpha_o12{2}(wc{end-Order:end-1});
        bpp1 = Coeffs.beta_o12{2}(wc{end-Order:end-1});
        gpp1 = Coeffs.gamma_o12{2}(wc{end-Order:end-1});
        %assert(gp(end) == gpp1(end));
    else % use family o26
        wflip = fliplr(w);
        if (Order == 1)
            apm1 = []; bpm1 = []; gpm1 = [];
        else
            apm1 = fliplr(Coeffs.alpha_o26{Order-1}(wflip(2:Order-1)));
            bpm1 = fliplr(Coeffs.beta_o26{Order-1}(wflip(2:Order-1)));
            gpm1 = fliplr(Coeffs.gamma_o26{Order-1}(wflip(2:Order-1)));
        end
        ap = fliplr(Coeffs.alpha_o26{Order}(wflip(2:Order)));
        bp = fliplr(Coeffs.beta_o26{Order}(wflip(2:Order)));
        gp = fliplr(Coeffs.gamma_o26{Order}(wflip(2:Order)));
        if (Order < Coeffs.maxOrder)
            app1 = fliplr(Coeffs.alpha_o26{Order+1}(wflip(2:Order+1)));
            bpp1 = fliplr(Coeffs.beta_o26{Order+1}(wflip(2:Order+1)));
            gpp1 = fliplr(Coeffs.gamma_o26{Order+1}(wflip(2:Order+1)));
            %assert(gp(end) == gpp1(end));
        else
            app1 = []; bpp1 = []; gpp1 = [];
        end
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
    
    gam = gp(end);
    gammas = zeros(k, 3);
    gammas(end-length(gpm1(1:end-1))+1:end,1) = gpm1(1:end-1)';
    gammas(end-length(gp(1:end-1))+1:end,2)   = gp(1:end-1)';
    gammas(end-length(gpp1(1:end-1))+1:end,3) = gpp1(1:end-1)';
end

function [Yp, YE, ELO, H, RejectStep, ISTATUS] = LIMMf_Onestep(Order, H, Y, F, k, IDX, Jacobian, dFdT, T, Coefficients, OdeFunction, OPTIONS, ISTATUS)

% a_(n+1) y_(n+1) + sum_i a_(n-i+1) y_(n-i+1) - h_n f(y_(n+1)) = 0

    [alpha, beta, gamma, gam] = Coefficients.getCoeffs(Coefficients, Order, H(IDX), k);
    
    methods = 1:3;
    if (Order == 1)
        methods = methods(2:end);
    end
    if (Order + 1 > nnz(H))
        methods = methods(1:end-1);
    end
    
    NVAR = length(Y(:,IDX(k)));
    Ynew = repmat(Y(:,IDX(k)), [1,3]);
    
    % build LHS matrices
    if (~OPTIONS.MatrixFree)
        LHS = speye(NVAR) - H(IDX(k)) * gam * Jacobian;
    else
        LHS = @(v)(v - H(IDX(k)) * gam * Jacobian(v));
    end
    
    %keyboard
    
    if (~OPTIONS.MatrixFree)
        RHS = -(Y(:,IDX) * alpha(:,methods)) + H(IDX(k)) * (F(:,IDX) * beta(:,methods)) + H(IDX(k)) * Jacobian * (Y(:,IDX) * gamma(:,methods));
        Ynew(:,methods) = LHS\RHS;
        ISTATUS.Nsol = ISTATUS.Nsol + 3;
    else
        for i = methods
            RHS = -(Y(:,IDX) * alpha(:,i)) + H(IDX(k)) * (F(:,IDX) * beta(:,i)) + H(IDX(k)) * Jacobian((Y(:,IDX) * gamma(:,i)));
            [dY, gmresFlag, ~, iter] = gmres(LHS, RHS, OPTIONS.GMRES_Restart, OPTIONS.GMRES_TOL, OPTIONS.GMRES_MaxIt, OPTIONS.GMRES_P);
            ISTATUS.Nsol = ISTATUS.Nsol + 1;

            if ( ~isempty(OPTIONS.GMRES_Restart) )
                vecCount = iter(2) + (OPTIONS.GMRES_Restart - 1)*iter(1);
            else
                vecCount = iter(2);
            end
            ISTATUS.Njac =  ISTATUS.Njac + vecCount;

            if( gmresFlag ~= 0 )
                resvec = abs(LHS(dY) - RHS);
                scalar = OPTIONS.AbsTol + OPTIONS.RelTol.*abs(RHS);
                if (norm(resvec./scalar) > sqrt(NVAR))
                    ISTATUS.Nsng = ISTATUS.Nsng + 1;
                    switch(gmresFlag)
                        case 1
                            warning('GMRES: iterated MAXIT times but did not converge');
                            break;
                        case 2
                            warning('GMRES: preconditioner M was ill-conditioned');
                            break;
                        case 3
                            warning('GMRES: stagnated (two consecutive iterates were the same)');
                            break;
                    end
                else
                    gmresFlag = 0;
                end
            end
            Ynew(:,i) = dY;
        end
    end
    
    Yp   = Ynew(:,2);
    
    % Compute error estimates
    YE = zeros(NVAR,3);
    ELO = ones(1,3);
    if ( nnz(H) < 2 ) % for the first step only
        assert(Order == 1, 'Order is not 1 for the first step.');
        [Ypp1, failureFlag, ISTATUS] = Rosenbrock2_onestep(H(IDX(k)), Y(:,IDX(k)), F(:,IDX(k)), Jacobian, T, OdeFunction, OPTIONS, ISTATUS);

        YE(:,1) = Yp - Ynew(:,1);
        YE(:,2) = Ypp1 - Yp;
        YE(:,3) = Yp - Ynew(:,3);
        
        ELO(1) = 1;
        ELO(2) = 2;
        ELO(3) = 1;
        
        if (failureFlag)
            warning('Rosenbrock2 failed to converge.');
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


