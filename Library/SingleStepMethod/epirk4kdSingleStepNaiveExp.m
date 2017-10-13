%% epirk4KSingleStepNaiveExp
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%
%% Input Parameters
%
%% Output Parameters
%
%% Description
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [y, yerr, ISTATUS] = epirk4kdSingleStepNaiveExp(y0, dt, rhsFun, jacFun, ...
                                    f, MatrixFree, NBasisVectors, ISTATUS, absTol, relTol, adaptiveKrylov)

    % Stages
    s = 3;

    % Coefficients computed using Mathematica

    % % Coefficients in front of Psi functions for all stages except last
    % a = [692665874901013/799821658665135, 0, 0;
    %     692665874901013/799821658665135, 3/4, 0];

    % % Final stage coefficients for 4th order method in front of Psi functions
    % b = [799821658665135/692665874901013, 352/729, 64/729];

    % % Final stage coefficients for 3rd order embedded method in front of Psi functions
    % b_hat = [799821658665135/692665874901013, 32/81, 0];

    % % Coefficients of matrices inside the exponentials
    % g = [3/4, 0, 0;
    %     3/4, 0, 0;
    %     1, 9/16, 9/16];

    % % Coefficients that take linear combinations of Phi functions
    % p = [692665874901013/799821658665135, 0 ,0;
    %     1,1,0;
    %     1,1,0];

    % Using the coefficients from the W-method to test its effectiveness as a K=method
    
    a = [ 0.228241829611716203969390857002, 0, 0;
         0.456483659223432407938781714004, 0.331616640633569500852166108254, 0];

    b = [1.00000000000000133922279642575, 2.0931591383832578214485343033, ...
         1.26239692579008044040805519461];

    b_hat = [1.00000000000000133922279642575, 2.0931591383832578214485343033, 1];


    g = [0, 0, 0;
         0.34706341174296320957640998939, 0.34706341174296320957640998939, ...
      0.34706341174296320957640998939;
        1.00000000000000000000000000000, 1.00000000000000000000000000000, ...
      1.00000000000000000000000000000];

    p = [0.99999999999999866077720357425, 0, 0;
         0, 2.09316041004385010039852804737, 0;
         0.99999999999999938123270102666, 0.99999999999999836439667546869, ...
         0.99999999999999794563711369222];
    
    % length of initial vector/autonomous system
    N = length(y0);

    % Stage vectors
    Y_s = zeros(N, s + 1);                                      % Y_i
    F_s = zeros(N, s + 1);                                      % F(Y_i)
    delta_s = zeros(N, s + 1);                                  % delta

    % Initialize Y_s with y0
    Y_s(:, 1) = y0;

    % Given s, construct a table of the combinations
    % j runs from 1 to s. The combinations however
    % upto j - 1 => s - 1
    C = zeros(s-1, s);
    for i = 1:s-1
        for j = 0:i
            C(i,j + 1) = nchoosek(i,j);
        end
    end

    % Given s, compute the p_tilde table
    p_tilde = zeros(s, 1);
    for j = 1:s
        for k = 1:j
            p_tilde(j) = p_tilde(j) + p(j, k)/factorial(k);
        end
    end

    % f is already computed outside. Hence commenting out next line.
    % f = rhsFun(Y_s(:,1));                   % Evaluate the rhsFun at y

    % J is computed inside ArnoldiAdapt. Hence commenting out next line.
    % J = jacFun(Y_s(:,1));                   % Evaluate the jacFun at y


    % Compute the Krylov basis matrices
    [V, H, M] = ArnoldiAdapt(jacFun, f, N, dt, MatrixFree, NBasisVectors);
    ISTATUS.Nkdim = ISTATUS.Nkdim + M^2;

    % Reduced Stage Vectors
    L_s = zeros(M, s + 1);                              % Lambda
    E_s = zeros(M, s + 1);                              % Eta

    F_s(:, 1) = f;                                      % F(Y_0)
    L_s(:, 1) = V' * Y_s(:, 1);                         % Lambda_0
    E_s(:, 1) = V' * F_s(:, 1);                         % Eta_0
    mu_0 = (Y_s(:, 1) - V * L_s(:, 1));                 % mu_0
    delta_s(:, 1) = (f - V * E_s(:, 1));                % delta_0

    % Stages run from 1 to s - 1.
    for i = 1 : s - 1
        L_s(:, i + 1) = L_s(:, 1);

        % Compute the Psi function only if a(i, 1) ~= 0
        if a(i, 1) ~= 0
            L_s(:, i + 1) = L_s(:, i + 1) + a(i, 1) * dt * ...
                PsiM(i, 1, H, M, dt, p_tilde(1), g, p) * E_s(:, 1);
        end

        for j = 2 : i
            if a(i, j) ~= 0
                d_j_1 = D_j_1(j - 1, H, M, E_s, L_s, C);        % Note that we can store and reuse previously
                L_s(:, i + 1) = L_s(:, i + 1) + ...             % computed d_j_1 for smaller j's. Not being done now.
                    a(i, j) * dt * ...
                    PsiM(i, j, H, M, dt, p_tilde(j), g, p) * d_j_1;
            end
        end

        Y_s(:, i + 1) = V * L_s(:, i + 1) + mu_0 + ...
            dt * a(i, 1) * p_tilde(1) * delta_s(:, 1);

        for j = 2 : i
            if a(i, j) ~= 0
                r_j_1 = R_j_1(j - 1, N, delta_s, C);
                Y_s(:, i + 1) = Y_s(:, i + 1) + ...
                    a(i, j) * dt * p_tilde(j) * r_j_1;
            end
        end

        F_s(:, i + 1) = rhsFun(Y_s(:, i + 1));                         % Evaluate the rhsFun at y_i
        E_s(:, i + 1) = V'*F_s(:, i + 1);                              % Eta_i
        delta_s(:, i + 1) = (F_s(:, i + 1) - V * E_s(:, i + 1));       % delta_0
    end


    % Final stage
    L_s(:, s + 1) = L_s(:, 1);

    % Initialize the error from the embedded method
    yerr_R = zeros(M, 1);                               % zeros in
                                                        % reduced
                                                        % subspace
    yerr = zeros(N, 1);                                 % zeros in
                                                        % full subspace


    % Compute the Psi function only if b(1) ~= 0
    psiTerm = dt * PsiM(s, 1, H, M, dt, p_tilde(1), g, p) * E_s(:, 1);
    L_s(:, s + 1) = L_s(:, s + 1) + b(1) * psiTerm;
    yerr_R = yerr_R + (b(1) - b_hat(1)) * psiTerm;

    for j = 2 : s
        d_j_1 = D_j_1(j - 1, H, M, E_s, L_s, C);        % Note that we can store and reuse previously
        psiTerm = dt *  PsiM(s, j, H, M, dt, p_tilde(j), g, p) * d_j_1;
        L_s(:, s + 1) = L_s(:, s + 1) + b(j) * psiTerm; % computed d_j_1 for smaller j's. Not being done now.
        yerr_R = yerr_R + (b(j) - b_hat(j)) * psiTerm;
    end

    psiTerm = dt * p_tilde(1) * delta_s(:, 1);
    Y_s(:, s + 1) = V * L_s(:, s + 1) + mu_0 + ...
             b(1) * psiTerm;
    yerr = V * yerr_R + (b(1) - b_hat(1)) * psiTerm;

    for j = 2 : s
        r_j_1 = R_j_1(j - 1, N, delta_s, C);
        psiTerm = dt * p_tilde(j) * r_j_1;
        Y_s(:, s + 1) = Y_s(:, s + 1) ...
            +  b(j) * psiTerm;
        yerr = yerr + (b(j) - b_hat(j)) * psiTerm;
    end

    % Get the solution at the next timestep
    y = Y_s(:, s + 1);
end


% This function computes d_(j - 1) for each j according to the definition
% (28) in the latex file.
% Need to pass j - 1 to the function instead of j.
% The function internally accounts for MATLAB indexing which starts at 1
% instead of 0
function d_j_1 = D_j_1(j_1, H, M, E_s, L_s, C)
    powers_of_negative_1 = 1;
    eta_sum = zeros(M, 1);
    lambda_sum = zeros(M, 1);
    for k = 0:j_1
        eta_sum = eta_sum + powers_of_negative_1 ...
            * C(j_1, k + 1) * E_s(:, j_1 - k + 1);
        lambda_sum = lambda_sum + powers_of_negative_1 ...
            * C(j_1, k + 1) * L_s(:, j_1 - k + 1);
        powers_of_negative_1 = powers_of_negative_1 * -1;
    end
    d_j_1 = eta_sum - H * lambda_sum;
end


% This function computes r_(j - 1) for each j according to the definition
% (29) in the latex file.
% Need to pass j - 1 to the function instead of j.
% The function internally accounts for MATLAB indexing which starts at 1
% instead of 0
function r_j_1 = R_j_1(j_1, N, delta_s, C)
    powers_of_negative_1 = 1;
    r_j_1 = zeros(N, 1);
    for k = 0:j_1
        r_j_1 = r_j_1 + powers_of_negative_1 ...
            * C(j_1, k + 1) * delta_s(:, j_1 - k + 1);
        powers_of_negative_1 = powers_of_negative_1 * -1;
    end
end

% In our formulation, the Psi function is dependent on
% i) parameter gij
% iii) the matrix H which stays the same.
function psi = PsiM(i, j, H, M, dt, p_tilde_j, g, p)
    % Note: that this function implements the simplified Psi defn
    % given in section 4 of EPIRK paper:
    % A new class of exponential propagation iterative methods of
    % Runge-Kutta type (EPIRK) - M. Tokman

    if g(i, j) ~= 0
        % Compute psi matrix in reduced space.
        psi = zeros(M);

        % intermediate phi_k for each k computed iteratively.
        % Using highPhi: renamed as phiKs
        [~, phi_ks] = phiKs(g(i,j) * dt * H, M, 1:j);

        for k = 1:j
           psi = psi + p(j, k) * phi_ks(:, :, k);
        end
    else
        % return the diagonal matrix with the p_tilde_j values along the
        % diagonal
        psi = diag(ones(M, 1) * p_tilde_j);
    end
end
