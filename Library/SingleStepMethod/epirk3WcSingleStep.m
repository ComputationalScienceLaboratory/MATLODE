%% epirk3WcSingleStepNaiveExp
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
function [y, yerr, ISTATUS] = epirk3WcSingleStep(y0, dt, rhsFun, jacFun, ...
                                    f, MatrixFree, NBasisVectors, ISTATUS, absTol, relTol, adaptiveKrylov)



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
    
    %----------------------------------------------------------------------
    % Set poor Tolerance level for Arnoldi - Ok for W method
    arnoldiTol = 1e-9;
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % Set the number of minimum basis vectors
    MBasisVectors = 1;
    %----------------------------------------------------------------------
    
    % f is already computed outside. Hence commenting out next line.
    % f = rhsFun(Y_s(:,1));                   % Evaluate the rhsFun at y

    % J is computed inside ArnoldiAdapt. Hence commenting out next line.
    % J = jacFun(Y_s(:,1));                   % Evaluate the jacFun at y

    % Stages run from 1 to s - 1. However, Y_s[1] = y_(i-1).
    % First stage
    Y_s(:, 2) = Y_s(:, 1);
    
    % Compute the Psi function only if a(i, 1) ~= 0
    % Term 2
    % Compute the Krylov basis matrices for column 1 (f)
    [V1, H1, M1] = ArnoldiAdapt(jacFun, f, N, max(g(:,1)) * dt, MatrixFree, ...
                             NBasisVectors, relTol, MBasisVectors, adaptiveKrylov);
    krySteps = krySteps + M1^2;
    
    [psi] = Psi(1, 1, f, N, dt, g, p, V1, H1, M1);
    Y_s(:, 2) = Y_s(:, 2) + a(1, 1) * dt * psi;

    % Second stage
    Y_s(:, 3) = Y_s(:, 1);

    % Compute the Psi function only if a(i, 1) ~= 0
    % Term 2
    [psi] = Psi(2, 1, f, N, dt, g, p, V1, H1, M1);
    Y_s(:, 3) = Y_s(:, 3) + a(2, 1) * dt * psi;
    % Term 3
    % Compute the Krylov basis matrices for column 2 (dr1)
    dr1 = dR(rhsFun, Y_s, 1, 1, jacFun, MatrixFree);
    [V2, H2, M2] = ArnoldiAdapt(jacFun, dr1, N, max(g(:,2)) * dt, MatrixFree, ...
                         NBasisVectors, relTol, MBasisVectors, adaptiveKrylov);
    krySteps = krySteps + M2^2;
    
    [psi] = Psi(2, 2, dr1, N, dt, g, p, V2, H2, M2);
    Y_s(:, 3) = Y_s(:, 3) + a(2, 2) * dt * psi;

    % Final stage
    Y_s(:, s + 1) = Y_s(:, 1);

    % Initialize the error vector computed using the embedded method.
    yerr = zeros(N, 1);

    % Compute the Psi function only if b(1) ~= 0
    % Term 2
    [psi] = Psi(s, 1, f, N, dt, g, p, V1, H1, M1);
    psiTerm = dt * psi;
    Y_s(:, s + 1) = Y_s(:, s + 1) + b(1) * psiTerm;
    yerr = yerr + (b(1) - b_hat(1)) * psiTerm;
    
    % Term 3
    [psi] = Psi(s, 2, dr1, N, dt, g, p, V2, H2, M2);
    psiTerm = dt * psi;
    Y_s(:, s + 1) = Y_s(:, s + 1) + b(2) * psiTerm;
    yerr = yerr + (b(2) - b_hat(2)) * psiTerm;
    
    % Term 4
    dr2 = dR(rhsFun, Y_s, 1, 2, jacFun, MatrixFree);
    % Compute the Krylov basis matrices for column 3 (dr2)
    [V3, H3, M3] = ArnoldiAdapt(jacFun, dr2, N, max(g(:,3)) * dt, MatrixFree, ...
                         NBasisVectors, relTol, MBasisVectors, adaptiveKrylov);
    krySteps = krySteps + M3^2;
    
    [psi] = Psi(s, 3, dr2, N, dt, g, p, V3, H3, M3);
    psiTerm = dt * psi;
    Y_s(:, s + 1) = Y_s(:, s + 1) + b(3) * psiTerm;
    yerr = yerr + (b(3) - b_hat(3)) * psiTerm;

    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/3;
    
    % Get the solution at the next timestep
    y = Y_s(:, s + 1);
end


% In our formulation, the Psi function is dependent on
% i) parameter gij
% ii) vector being multiplied by Psi function
% iii) the matrix A0 which stays the same.
function [psi] = Psi(i, j, v, N, dt, g, p, V, H, M)
    % Note: that this function implements the simplified Psi defn
    % given in section 4 of EPIRK paper:
    % A new class of exponential propagation iterative methods of
    % Runge-Kutta type (EPIRK) - M. Tokman

    % Refer to the standard implementation (section 5.1) in
    % Exponential-Krylov methods for ODEs - Tranquilli, Sandu
    % (doi:10.1016/j.jcp.2014.08.013)
    % paper for details on why unit vector e1 appears.
    % Refer to Theorem 1 in ExpoKit on how HBar is constructed
    % Note: 1) that Tau = 1 so as to avoid all the Tau^2 Tau^3 etc.
    % The multiplier gij has already multiplied HBar.
    % 2) dt between the Psi and vector is not actually multiplied,
    % do it outside the method or here. Outside better for clarity.
    if g(i, j) ~= 0
        psi = zeros(N, 1);
        e1 = [1; zeros(M - 1, 1)];
        Hbar = [g(i,j) * dt * H e1 zeros(M, j-1);zeros(j-1, M+1) eye(j - 1);zeros(1, M+j)];
        expHbar = expm(Hbar);                    % Tau = 1
        normv = norm(v, 2);

        for k = 1:j
           psi = psi + p(j, k) *  normv * V * expHbar(1:M, M + k);
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
        r = rhsFun(y) - rhsFun(y0) - A0 * (y - y0);         % Need to find if
                                                            % A0 has to be
                                                            % the actual
                                                            % Jacobian or
                                                            % Krylov subspace
                                                            % approximation.
    else
        r = rhsFun(y) - rhsFun(y0) - A0(y - y0);            % Need to find if
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
