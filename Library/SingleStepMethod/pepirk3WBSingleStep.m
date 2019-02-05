function [y, yerr, ISTATUS] = pepirk3WBSingleStep(y0, dt, rhsFun1, rhsFun2, ...
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Stages
    s = [3,3];
    
    % Partitions
    ps = 2;
    
    % Coefficients for third-order EPIRK method
    A11 = [[-2/3, 0, 0]; [-2/3, 643210555/321605294, 0]; [0, 0, 0]];
    A12 = [[-77257927/231736567, 0, 0];  [-77257927/231736567, 221995883/556342956, 0]; [0, 0, 0]];  
    A21 = [[46177625/138532819, 0, 0]; [46177625/138532819, -276294928/448265855, 1]; [0, 0, 0]];
    A22 = [[2/3, 0, 0]; [2/3, -2097151/1048576, 1]; [0, 0, 0]];
    A = {A11, A12; A21, A22};
    
    
    % Coupling happens through the A coefficients
    G11 = [[15905698/118882405, 0, 0]; [47844971/181848383, 0, 0]; [0, 263476/712854121, 5873/47718781]];   
    G12 = [[7866/337976873, 0, 0]; [11957/261254373, 0, 0]; [0.1, 0.1, 0.1]];
    G21 = [[151156669/742595986, 0, 0]; [111733984/275857671, 0, 0]; [0.1, 0.1, 0.1]];
    G22 = [[0, 0, 0]; [0, 0, 0]; [0, 57491/77959816, 2442/336507713]];
    G = {G11, G12; G21, G22};
    
    % Coupling happens through the A coefficients
    P11 = [[-1, 0, 0]; [-90974766/422643103, 15713557/145038387, 0]; [139909573/248730863, -355468325/226721754, -279447475/279556943]];   
    P12 = [[-534076380/267081073, 0, 0]; [-306646172/379808307, -1391/610013106, 0]; [1, 1, 1]];
    P21 = [[291908440/145954279, 0, 0]; [341330517/637547050, -10114/98573483, 0]; [1, 1, 1]];
    P22 = [[1, 0, 0]; [21675617/200502754, 27071061/237997009, 0]; [429818508/221501921, 813095776/411099863, -541652055/542849144]];
    P = {P11, P12; P21, P22};
        
    B1 =  [-1, 36471146/225920009, 92167149/46087892];
    B2 =  [1, -15699629/337585791, -69127796/252097195];
    B = {B1,B2};
    
    BH1 =  [-1, 1, 49043754/20888203];
    BH2 =  [1, 1, -25364743/119799726];
    BH  = {BH1,BH2};
    
    % Size of problem
    N = length(y0);
    
    % Stage vectors
    Y_s = zeros(N, max(s(:)), ps);
    F_s = zeros(N, 1, ps);          % Only the function value corresponding to Y_0.
    
    % optional arguments
    % Set the minimum number of basis vectors for the method
    if(~exist('MBasisVectors','var'))
        MBasisVectors = 1;
    end
    
    % create an array of function/matrix handles
    rhsFuns = {rhsFun1, rhsFun2};
    jacFuns = {jac1, jac2};


    % krylov steps
    krySteps = 0;
    
    % psi_count
    psi_count = 0;
    
    % Jacobians or their matrix-vector product functions
    % for all partitions.
    JacJ = cell(ps,1);
    
    % current time
    Y_s(:, 1, 1) = y0;   % partition 1
    Y_s(:, 1, 2) = y0;   % partition 2
    F_s(:, 1, 1) =  f1_0; % partition 1  % reuse if sent in [else rhsFun1(Y_s(:, 1, 1));]
    F_s(:, 1, 2) =  f2_0; % partition 2  % reuse if sent in [else rhsFun2(Y_s(:, 1, 2));]
    
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
    for i = 1:max(s(:)) - 1                             % Stage Index upto last but one
        % Do for each partition
        % First Exponential
        % TODO: This can be refactored to accept an array of
        % function handles and use that instead
        for j = 1:ps                                    % Partition Index
            % if the current stage index is less
            % than the total number of stages in the
            % j'th partition then do the following
            if i <= s(j) - 1
                Y_s(:, i + 1, j) = y0;
                
                % computes all the mixed first term
                for k = 1:ps
                    %term2 = A{j,k}(i, 1) * dt * PsiExact(i, 1, JacJ{k}, F_s(:, 1, k),...
                    %    N, dt, G{j,k}, P{j,k}, MatrixFree, NBasisVectors, krySteps, relTol);
                    term = A{j,k}(i, 1) * dt * Psi(i, 1, JacJ{k}, F_s(:, 1, k),...
                        N, dt, G{j,k}, P{j,k}, MatrixFree, NBasisVectors, krySteps, relTol);
                    psi_count = psi_count + 1;
                    Y_s(:, i + 1 , j) = Y_s(:, i + 1, j) + term;
                end
                
                for k = 1:ps
                    for l = 2:min(i, s(j) - 1)
                        dr = dR(rhsFuns{k}, Y_s(:, :, k), 1, l - 1, JacJ{k}, MatrixFree);
                        %term2 = A{j,k}(i, l) * dt * PsiExact(i, l, JacJ{k}, dr,...
                        %    N, dt, G{j,k}, P{j,k}, MatrixFree, NBasisVectors, krySteps, relTol);
                        term = A{j,k}(i, l) * dt * Psi(i, l, JacJ{k}, dr,...
                            N, dt, G{j,k}, P{j,k}, MatrixFree, NBasisVectors, krySteps, relTol);
                        psi_count = psi_count + 1;
                        Y_s(:, i + 1, j) = Y_s(:, i + 1, j) + term;
                    end
                end
            end
        end
    end
    
    y = y0;
    yerr = zeros(size(y0));
    % computes all the mixed first term
    for k = 1:ps
        %term2 = B{k}(1) * dt * PsiExact(s(k), 1, JacJ{k}, F_s(:, 1, k), ...
        %    N, dt, G{k,k}, P{k,k}, MatrixFree, NBasisVectors, krySteps, relTol);
        term = B{k}(1) * dt * Psi(s(k), 1, JacJ{k}, F_s(:, 1, k), ...
            N, dt, G{k,k}, P{k,k}, MatrixFree, NBasisVectors, krySteps, relTol);
        psi_count = psi_count + 1;
        y = y + term;
        yerr = yerr + term/ B{k}(1) * (B{k}(1) - BH{k}(1));
    end
    
    for k = 1:ps
        for l = 2:s(k)
            dr = dR(rhsFuns{k}, Y_s(:, :, k), 1, l - 1, JacJ{k}, MatrixFree);
            %term2 = B{k}(l) * dt * PsiExact(s(k), l, JacJ{k}, dr, ...
            %    N, dt, G{k,k}, P{k,k}, MatrixFree, NBasisVectors, krySteps, relTol);
            term = B{k}(l) * dt * Psi(s(k), l, JacJ{k}, dr, ...
                N, dt, G{k,k}, P{k,k}, MatrixFree, NBasisVectors, krySteps, relTol);
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
    if MatrixFree
        r = rhsFun(y) - rhsFun(y0); %- jacFun(y - y0);       % Need to find if
    else
        r = rhsFun(y) - rhsFun(y0);% - jacFun * ( y - y0);   % Need to find if
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
