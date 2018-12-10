%% epirk3WDSingleStepNaiveExp
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
function [y, yerr, ISTATUS] = epirk3WcTriDSingleStepNaiveExp(y0, dt, rhsFun, jacFun, ...
                                    f, MatrixFree, NBasisVectors, ISTATUS, absTol, relTol, adaptiveKrylov, symmjac, MBasisVectors)

    disp('here second');

    % Stages
    s = 3;

    % Coefficients computed using Mathematica
    a = [ 282/311,     0, 0;
          294/311, -7/94, 0];

    b = [1, -3421/987, -622/105];

    b_hat = [1, 13/9, 1];


    g = [1/5, 0,   0;
         1/8, 1/8, 0;
         1,   1,   1];

    p = [1,     0, 0;
         1/2, 1/2, 0;
         1/3, 1/3, 1/3];

    % length of initial vector/autonomous system
    N = length(y0);

    % Stage vectors
    Y_s = zeros(N, s + 1);                                      % Y_i


    % Initialize Y_s with y0
    Y_s(:, 1) = y0;                                             % Stage vectors

    krySteps = 0;
    psiCount = 0;
   
    %----------------------------------------------------------------------
    % Set poor Tolerance level for Arnoldi
    arnoldiTol = 1e-12;
    %----------------------------------------------------------------------

    % f is already computed outside. Hence commenting out next line.
    % f = rhsFun(Y_s(:,1));                   % Evaluate the rhsFun at y

    % J is computed inside ArnoldiAdapt. Hence commenting out next line.
    % J = jacFun(Y_s(:,1));                   % Evaluate the jacFun at y

    % Stages run from 1 to s - 1. However, Y_s[1] = y_(i-1).
    for i = 1 : s - 1
        Y_s(:, i + 1) = Y_s(:, 1);

        % Compute the Psi function only if a(i, 1) ~= 0
        if a(i, 1) ~= 0
            [psi, krySteps] = Psi(i, 1, jacFun, f, N, ...
                                  dt, g, p, MatrixFree, NBasisVectors, krySteps, arnoldiTol);
            psiCount = psiCount + 1;
            Y_s(:, i + 1) = Y_s(:, i + 1) + a(i, 1) * dt * psi;
        end

        for j = 2: i
            if a(i, j) ~= 0
                dr = dR(rhsFun, Y_s, 1, j - 1, jacFun, MatrixFree);
                [psi, krySteps] = Psi(i, j, jacFun, dr, N, dt,  g, p, ...
                                      MatrixFree, NBasisVectors, krySteps, arnoldiTol);
                psiCount = psiCount + 1;
                Y_s(:, i + 1) = Y_s(:, i + 1) + a(i, j) * dt * psi;
            end
        end
    end

    % Final stage
    Y_s(:, s + 1) = Y_s(:, 1);

    % Initialize the error vector computed using the embedded method.
    yerr = zeros(N, 1);

    % Compute the Psi function only if b(1) ~= 0
    if b(1) ~= 0
        [psi, krySteps] = Psi(s, 1, jacFun, f, N, dt, g, p, MatrixFree,...
                              NBasisVectors, krySteps, arnoldiTol);
        psiCount = psiCount + 1;
        psiTerm = dt * psi;
        Y_s(:, s + 1) = Y_s(:, s + 1) + b(1) * psiTerm;
        yerr = yerr + (b(1) - b_hat(1)) * psiTerm;
    end

    for j = 2: s
        if b(j) ~= 0
            dr = dR(rhsFun, Y_s, 1, j - 1, jacFun, MatrixFree);
            [psi, krySteps] = Psi(s, j, jacFun, dr, N, dt, g, p, MatrixFree, ...
                                  NBasisVectors, krySteps, arnoldiTol);
            psiCount = psiCount + 1;
            psiTerm = dt * psi;
            Y_s(:, s + 1) = Y_s(:, s + 1) + b(j) * psiTerm;
            yerr = yerr + (b(j) - b_hat(j)) * psiTerm;
        end
    end

    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/psiCount;

    % Get the solution at the next timestep
    y = Y_s(:, s + 1);
end


% In our formulation, the Psi function is dependent on
% i) parameter gij
% ii) vector being multiplied by Psi function
% iii) the matrix A0 which stays the same.
function [psi, krySteps] = Psi(i, j, A0, v, N, dt, g, p, MatrixFree, NBasisVectors, krySteps, Tol)
    % Note: that this function implements the simplified Psi defn
    % given in section 4 of EPIRK paper:
    % A new class of exponential propagation iterative methods of
    % Runge-Kutta type (EPIRK) - M. Tokman

    if(~exist('Tol','var'))
        arnoldiTol = 1e-12;
    else
        arnoldiTol = Tol;
    end
    
    if g(i, j) ~= 0
        psi = zeros(N, 1);
        
        if (~MatrixFree)
            d = diag(A0);
        else
            N = size(A0, 1);
            I = eye(N);
            D = zeros(N);
            for i = 1:N
                D(:, i) = A0(I(:,i));
            end
            d = diag(D);
        end
        
        psiM = g(i,j)*dt*d;
        psiM = exp(psiM);                                           % phi_0

        for k = 1:j
           psiM = (psiM - 1/factorial(k - 1))./(g(i,j) * dt * d);   % phi_k
           psi = psi + p(j, k) *  psiM .* v;
        end
    else
        p_tilde = 0;
        for k = 1:j
           p_tilde  = p_tilde + p(j, k)/factorial(k);
        end
        psi = p_tilde * v;
    end
end


% R function.
function r = R(rhsFun, y, y0, A0, MatrixFree)
    if (~MatrixFree)
        r = rhsFun(y) - rhsFun(y0) - diag(A0) .* (y - y0);  % Need to find if
                                                            % A0 has to be
                                                            % the actual
                                                            % Jacobian or
                                                            % Krylov subspace
                                                            % approximation.
    else
        N = size(A0, 1);
        I = eye(N);
        D = zeros(N);
        for i = 1:N
            D(:, i) = A0(I(:,i));
        end
        
        r = rhsFun(y) - rhsFun(y0) - diag(D) .* (y - y0);   % Need to find if
                                                            % A0 has to be
                                                            % the actual
                                                            % Jacobian or
                                                            % Krylov subspace
    end
end

% Forward Difference
% Note: It is likely that recompute and dr are not being utilized.
% Matlab copies on write. Do the changes then reflect in the original?
function dr = dR(rhsFun, Y_s, i, j, A0, MatrixFree)
    if j == 1
            dr = R(rhsFun, Y_s(:,i+1), Y_s(:,1), A0, MatrixFree) - ...
                 R(rhsFun, Y_s(:,i), Y_s(:,1), A0, MatrixFree);
    elseif j > 1
        dr = dR(rhsFun, Y_s, i+1, j-1, A0, MatrixFree) ...
            - dR(rhsFun, Y_s, i, j-1, A0, MatrixFree);
    else
        error('Incorrect j value');
    end
end
