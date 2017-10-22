%% epirk4kbSingleStep
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
function [y, yerr, ISTATUS] = epirk4kbSingleStep(y0, dt, rhsFun, jacFun, ...
                                    f, MatrixFree, NBasisVectors, ISTATUS, absTol, relTol, adaptiveKrylov, symmjac);
    % Stages
    s = 3;

    % Coefficients computed using Mathematica

    % Coefficients in front of Psi functions for all stages except last
    a = [1, 0, 0;
         1, 1, 0];

    % Final stage coefficients for 4th order method in front of Psi functions
    b = [4/3, 112/243, 1];

    % Final stage coefficients for 3rd order embedded method in front of Psi functions
    b_hat = [4/3, 80/243, -1];

    % Coefficients of matrices inside the exponentials
    g = [3/4,  0,   0;
        3/4,   3/4, 0;
        1,     3/4, 3/4];

    % Coefficients that take linear combinations of Phi functions
    p = [3/4, 0 , 0;
        1,    1,  0;
        1,    -962/243, 524/81];

    % length of initial vector/autonomous system
    N = length(y0);

    % Stage vectors
    Y_s = zeros(N, s + 1);                                  % Y_i
    F_s = zeros(N, s + 1);                                  % F(Y_i)
    delta_s = zeros(N, s + 1);                              % delta
    yerr = zeros(N, 1);                                     % error

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

    % Minimum number of basis vectors needed for 4th order convergence
    MBasisVectors = 4;

    % Set a high arnoldi tolerance
    arnoldiTol = 1e-12;
    
    % Compute the Krylov basis matrices
    [V, H, M] = ArnoldiAdapt(jacFun, f, N, dt, MatrixFree, NBasisVectors, relTol, MBasisVectors, adaptiveKrylov, symmjac);
    
    ISTATUS.Nkdim = ISTATUS.Nkdim + M^2;

    % Compute the LU factorization of H
    [LH, UH, pH] = lu(H,'vector');

    % Need phi1(3/4), phi2(3/4), phi3(3/4)
    [~, phi_ks] = phiKs(3/4 * dt * H, M, 1:3, LH, UH, pH, dt, 3/4);
    phi1_3_4 = phi_ks(:, :, 1);
    phi2_3_4 = phi_ks(:, :, 2);
    phi3_3_4 = phi_ks(:, :, 3);
    % Need phi1(1)
    [~, phi_ks] = phiKs(dt * H, M, 1, LH, UH, pH, dt, 1);
    phi1_1 = phi_ks(:, :, 1);

    % Reduced Stage Vectors
    L_s = zeros(M, s + 1);                              % Lambda
    E_s = zeros(M, s + 1);                              % Eta

    % Initialize the error from the embedded method
    yerr_R = zeros(M, 1);                               % zeros in reduced subspace


    F_s(:, 1) = f;                                      % F(Y_0)
    L_s(:, 1) = V' * Y_s(:, 1);                         % Lambda_0
    E_s(:, 1) = V' * F_s(:, 1);                         % Eta_0
    mu_0 = (Y_s(:, 1) - V * L_s(:, 1));                 % mu_0
    delta_s(:, 1) = (f - V * E_s(:, 1));                % delta_0

    % i = 1, Stage 1
    L_s(:, 2) = L_s(:, 1);
    L_s(:, 2) = L_s(:, 2) + a(1, 1) * dt * p(1,1) * phi1_3_4 * E_s(:, 1);

    % Compute the full-stage vector for first stage

    Y_s(:, 2) = V * L_s(:, 2) + mu_0 + a(1, 1) * dt * p_tilde(1) * delta_s(:, 1);

    F_s(:, 2) = rhsFun(Y_s(:, 2));
    E_s(:, 2) = V'*F_s(:, 2);
    delta_s(:, 2) = (F_s(:, 2) - V * E_s(:, 2));

    % i = 2, Stage = 2
    L_s(:, 3) = L_s(:, 1);
    L_s(:, 3) = L_s(:, 3) + a(2, 1) * dt * p(1,1) * phi1_3_4 * E_s(:, 1);

    d_2_1 = D_j_1(1, H, M, E_s, L_s, C);

    L_s(:, 3) = L_s(:, 3) + a(2, 2) * dt * (p(2,1) * phi1_3_4 + p(2,2) * phi2_3_4) * d_2_1;

    r_2_1 = R_j_1(1, N, delta_s, C);

    % Compute the full-stage vector for second stage
    Y_s(:, 3) = V * L_s(:, 3) + mu_0 + a(2, 1) * dt * p_tilde(1) * delta_s(:, 1) ...
                + a(2, 2) * dt * p_tilde(2) * r_2_1;
    F_s(:, 3) = rhsFun(Y_s(:, 3));
    E_s(:, 3) = V'*F_s(:, 3);
    delta_s(:, 3) = (F_s(:, 3) - V * E_s(:, 3));

    % i = 3, Stage = 3
    L_s(:, 4) = L_s(:, 1);

    psiTerm = dt * p(1,1) * phi1_1 * E_s(:, 1);
    L_s(:, 4) = L_s(:, 4) + b(1) * psiTerm;
    yerr_R = yerr_R + (b(1) - b_hat(1)) * psiTerm;

    psiTerm = dt * (p(2, 1) * phi1_3_4 + p(2, 2) * phi2_3_4) * d_2_1;
    L_s(:, 4) = L_s(:, 4) + b(2) * psiTerm;
    yerr_R = yerr_R + (b(2) - b_hat(2)) * psiTerm;

    d_3_1 = D_j_1(2, H, M, E_s, L_s, C);

    psiTerm = dt * (p(3, 1) * phi1_3_4 + p(3, 2) * phi2_3_4 + p(3, 3) * phi3_3_4)  * d_3_1;
    L_s(:, 4) = L_s(:, 4) + b(3) * psiTerm;
    yerr_R = yerr_R + (b(3) - b_hat(3)) * psiTerm;

    % Compute the full-stage vector for last stage
    psiTerm = dt * p_tilde(1) * delta_s(:, 1);
    Y_s(:, 4) = V * L_s(:, 4) + mu_0 + b(1) * psiTerm;
    yerr = V * yerr_R + (b(1) - b_hat(1)) * psiTerm;

    psiTerm = dt * p_tilde(2) * r_2_1;
    Y_s(:, 4) = Y_s(:, 4) + b(2) * psiTerm;
    yerr = yerr + (b(2) - b_hat(2)) * psiTerm;

    r_3_1 = R_j_1(2, N, delta_s, C);

    psiTerm = dt * p_tilde(3) * r_3_1;
    Y_s(:, 4) = Y_s(:, 4) + b(3) * psiTerm;
    yerr = yerr + (b(3) - b_hat(3)) * psiTerm;

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
