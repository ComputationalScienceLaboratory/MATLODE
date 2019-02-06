function [y, yerr, ISTATUS] = psepirk3bSingleStep(y0, dt, rhsFun1, rhsFun2, ...
        jac1, jac2, f1_0, f2_0, MatrixFree, NBasisVectors, ISTATUS, ...
        absTol, relTol, adaptiveKrylov, symmjac, MBasisVectors)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % EPIRK-like framework based GARK type method.
    %
    %
    % Output
    % y:     Solution at final time
    % yerr:  Error vector
    %
    % Input
    % y0     : Initial condition
    % dt     : timestep size
    % Tspan  : time interval over which solution is desired
    % rhsFun1: F^1, the first piece of the right hand side (first stiff piece to be treated exponentially)
    % rhsFun2: F^2, the second piece of the right hand side (second stiff piece to be treated exponentially)
    % jacFun1: Linear piece of the first function to be used in the exponential term
    % jacFun2: Linear piece of the second function to be used in the exponential term
    %
    % Method implemented for autonomous system
    % Mathematica file: EPIRK_Method15.nb
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Stages
    s = 3;
    
    % Partitions
    ps = 2;
    
    A1 = [[-406771190822767/1269704574921366, 0, 0]; ...
        [1141343580746374/2770232365403325, -270901875227597/1180295506525702, 0]; ...
        [0, 0, 0]];
    
    A2 = [[-406771190822767/1269704574921366, 0, 0]; ...
        [1141343580746374/2770232365403325, -270901875227597/1180295506525702, 0]; ...
        [0, 0, 0]];
    
    B1 = [791160475717460/1280186143788623, 411329415142993/512163583980334, 1410289361848454/2446680005682115];
    B2 = [791160475717460/1280186143788623, 411329415142993/512163583980334, 1410289361848454/2446680005682115];
    
    G1 = [[128087174255567/419295332948579, 0, 0]; ...
        [2/3, 1/3, 0]; ...
        [1, 1/10, 360563110168627/1057763542537753]];
    
    G2 = [[128087174255567/419295332948579, 0, 0]; ...
        [2/3, 1/3, 0]; ...
        [1, 1/10, 360563110168627/1057763542537753]];
    
    P1 = [[1280186143788623/1582320951434920, 0, 0]; ...
        [542606543994981/727846286011843, 334677002077769/888394725117090, 0]; ...
        [1537951833313021/2936537612285079, 174454083236061/758739097109116, 93287345944837/1308035706514201]];
    
    P2 = [[1280186143788623/1582320951434920, 0, 0]; ...
        [542606543994981/727846286011843, 334677002077769/888394725117090, 0]; ...
        [1537951833313021/2936537612285079, 174454083236061/758739097109116, 93287345944837/1308035706514201]];
    
    BH1 = [427455192538559/562171059204137, 382776428251149/390415271373872, 376136126102442/532899712497139];
    BH2 = [318902647845412/670465682883183, 1144042149029406/1168515417076291, 382460571524017/636059012992202];
    
    % Coefficients for third-order EPIRK method
    A = {A1, A2};
    
    
    G = {G1, G2};
    
    
    P = {P1, P2};
    
    % Final Stage
    B = {B1, B2};
    
    % Embedded Stage
    BH = {BH1, BH2};
    
    % Size of problem
    N = length(y0);
    
    % Stage vectors
    Y_s = zeros(N, max(s(:)));
    
    % optional arguments
    % Set the minimum number of basis vectors for the method
    if(~exist('MBasisVectors','var'))
        MBasisVectors = 1;
    end
    
    % create an array of function/matrix handles
    jacFuns = {jac1, jac2};
    
    % create a function handle for full F
    F = @(y) (rhsFun1(y) + rhsFun2(y));
    
    % krylov steps
    krySteps = 0;
    
    % psi_count
    psi_count = 0;
    
    % Jacobians or their matrix-vector product functions
    % for all partitions.
    JacJ = cell(ps,1);
    
    % current time
    Y_s(:, 1)    =  y0;
    F_0          =  f1_0 + f2_0;
    
    for j = 1:ps
        if (~MatrixFree)
            JacJ{j} = jacFuns{j};
        else
            % A partial jacobian function with the time
            % and state vector embedded is created in parent script
            JacJ{j} = @(v) jacFuns{j}(v);
        end
    end
    
    % Do for each stage
    for i = 1:s - 1                             % Stage Index upto last but one
        Y_s(:, i + 1) = y0;
        
        for j = 1:ps
            term = A{j}(i, 1) * dt * Psi(i, 1, JacJ{j}, F_0,...
                N, dt, G{j}, P{j}, MatrixFree, NBasisVectors, krySteps, relTol);
            Y_s(:, i + 1) = Y_s(:, i + 1) + term;
        end
        
        for j = 1:ps
            for l = 2:min(i, s - 1)
                dr = dR(F, Y_s, 1, l - 1, JacJ{j}, MatrixFree);
                term = A{j}(i, l) * dt * Psi(i, l, JacJ{j}, dr,...
                    N, dt, G{j}, P{j}, MatrixFree, NBasisVectors, krySteps, relTol);
                Y_s(:, i + 1) = Y_s(:, i + 1) + term;
            end
        end
    end
    
    y = y0;
    yerr = zeros(size(y0));
    % computes all the mixed first term
    for k = 1:ps
        %term2 = B{k}(1) * dt * PsiExact(s, 1, JacJ{k}, F_0, ...
        %    N, dt, G{k}, P{k}, MatrixFree, NBasisVectors, krySteps, relTol);
        term = B{k}(1) * dt * Psi(s, 1, JacJ{k}, F_0, ...
            N, dt, G{k}, P{k}, MatrixFree, NBasisVectors, krySteps, relTol);
        psi_count = psi_count + 1;
        y = y + term;
        yerr = yerr + term/ B{k}(1) * (B{k}(1) - BH{k}(1)); % Remove divide by bk1
    end
    
    for k = 1:ps
        for l = 2:s
            dr = dR(F, Y_s, 1, l - 1, JacJ{k}, MatrixFree);
            %term2 = B{k}(l) * dt * PsiExact(s, l, JacJ{k}, dr, ...
            %    N, dt, G{k}, P{k}, MatrixFree, NBasisVectors, krySteps, relTol);
            term = B{k}(l) * dt * Psi(s, l, JacJ{k}, dr, ...
                N, dt, G{k}, P{k}, MatrixFree, NBasisVectors, krySteps, relTol);
            psi_count = psi_count + 1;
            y = y + term;
            yerr = yerr + term/ B{k}(l) * (B{k}(l) - BH{k}(l));
        end
    end
    
    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/psi_count;
end


% In our formulation, the Psi function is dependent on
% i) parameter gij
% ii) vector being multiplied by Psi function
% iii) the matrix A0 which stays the same.
function [phi1] = Phi_1(gij, v, dt, V, H, M)
    % Refer to the standard implementation (section 5.1) in
    % Exponential-Krylov methods for ODEs - Tranquilli, Sandu
    % (doi:10.1016/j.jcp.2014.08.013)
    % paper for details on why unit vector e1 appears.
    % Refer to Theorem 1 in ExpoKit on how HBar is constructed
    % Note: 1) that Tau = 1 so as to avoid all the Tau^2 Tau^3 etc.
    % The multiplier gij has already multiplied HBar.
    % 2) dt between the Psi and vector is not actually multiplied,
    % do it outside the method or here. Outside better for clarity.
    if gij ~= 0
        e1 = [1; zeros(M - 1, 1)];
        Hbar = [gij * dt * H e1;zeros(1, M+1)];
        expHbar = expm(Hbar);                    % Tau = 1
        normv = norm(v, 2);
        phi1 = normv * V * expHbar(1:M, M + 1);
    else
        phi1 = v;
    end
end

function [psiV] = PsiExact(i, j, A0, v, N, dt, g, p, MatrixFree, NBasisVectors, krySteps, Tol)
    % psi_ij(z) = psi_j(z) = sum(pjk * phi_k(z)) for k = 1 .. j
    % compute all the necessary phi_ks
    
    if (MatrixFree)
        AA = A0(eye(N));
    else
        AA = A0;
    end
    
    [~, phi_ks] = phiKs(AA * g(i, j) * dt, N, 1:j);
    
    psiV = zeros(N, 1);
    
    for k = 1:j
        psiV = psiV + p(j, k) * phi_ks(:, :, j) * v;
    end
end

% In our formulation, the Psi function is dependent on
% i) parameter gij
% ii) vector being multiplied by Psi function
% iii) the matrix A0 which stays the same.
function [psiV, krySteps] = Psi(i, j, A0, v, N, dt, g, p, MatrixFree, NBasisVectors, krySteps, Tol)
    % Note: that this function implements the simplified Psi defn
    % given in section 4 of EPIRK paper:
    % A new class of exponential propagation iterative methods of
    % Runge-Kutta type (EPIRK) - M. Tokman
    
    if(~exist('Tol','var'))
        arnoldiTol = 1e-12;
    else
        arnoldiTol = Tol;
    end
    
    % if v is zero, return 0
    if(norm(v) == 0)
        psiV = v;
        return
    end
    
    % Refer to the standard implementation (section 5.1) in
    % Exponential-Krylov methods for ODEs - Tranquilli, Sandu
    % (doi:10.1016/j.jcp.2014.08.013)
    % paper for details on why unit vector e1 appears.
    % Refer to Theorem 1 in ExpoKit on how HBar is constructed
    % Note: 1) that Tau = 1 so as to avoid all the Tau^2 Tau^3 etc.
    % The multiplier gij has already multiplied HBar.
    % 2) dt between the Psi and vector is not actually multiplied,
    % do it outside the method or here. Outside better for clarity.
    % TODO: Check norm of A?
    if g(i, j) ~= 0
        psiV = zeros(N, 1);
        
        % Compute the Krylov basis matrices
        [V, H, M] = ArnoldiAdapt(A0, v, N, g(i, j) * dt, MatrixFree, ...
            NBasisVectors, arnoldiTol);
        
        krySteps = krySteps + M^2;
        e1 = [1; zeros(M - 1, 1)];
        Hbar = [g(i,j) * dt * H e1 zeros(M, j-1);zeros(j-1, M+1) eye(j - 1);zeros(1, M+j)];
        expHbar = expm(Hbar);                    % Tau = 1
        normv = norm(v, 2);
        
        for k = 1:j
            psiV = psiV + p(j, k) *  normv * V * expHbar(1:M, M + k);
        end
    else
        p_tilde = 0;
        for k = 1:j
            p_tilde  = p_tilde + p(j, k)/factorial(k);
        end
        psiV = p_tilde * v;
    end
end


% R function.
function r = R(rhsFun, y, y0, jacFun, MatrixFree)
    if ~MatrixFree
        r = rhsFun(y) - rhsFun(y0) - jacFun * (y - y0);
    else
        r = rhsFun(y) - rhsFun(y0) - jacFun(y - y0);
    end
end

% Forward Difference
% Note: It is likely that recompute and dr are not being utilized.
% Matlab copies on write. Do the changes then reflect in the original?
function dr = dR(rhsFun, Y_s, i, j, jacFun, MatrixFree)
    if j == 1
        dr = R(rhsFun, Y_s(:,i+1), Y_s(:,1), jacFun, MatrixFree) - ...
            R(rhsFun, Y_s(:,i), Y_s(:,1), jacFun, MatrixFree);
    elseif j > 1
        dr = dR(rhsFun, Y_s, i+1, j-1, jacFun, MatrixFree) ...
            - dR(rhsFun, Y_s, i, j-1, jacFun, MatrixFree);
    else
        error('Incorrect j value');
    end
end
