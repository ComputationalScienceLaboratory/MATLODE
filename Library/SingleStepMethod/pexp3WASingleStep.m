function [y, yerr, ISTATUS] = pexp3WASingleStep(y0, dt, rhsFun1, rhsFun2, ...
        jacFun1, jacFun2, f1_0, f2_0, MatrixFree, NBasisVectors, ISTATUS, ...
        absTol, relTol, adaptiveKrylov, symmjac, NReactants, Autonomous, MBasisVectors)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % EXP framework based GARK type method.
    %
    %
    % Output
    % Y: Solution at all timesteps
    % y: Solution at final time
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
% Mathematica file: expK_framework_embedded_method_nosubsuperscripts_withdefaults.nb
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Stages
    s = [3,4];
    
    % Partitions
    ps = 2;
    
    % Coefficients for third-order ExpEx method
    % There is a discrepancy between the A's in the MATLAB file below and the one in Mathematica (depending on how code is written)
    % check with the formulation and adjust as necessary
    
   
    A11 = [[0, 0, 0];...
           [1, 0, 0]; ...
           [1/9, 2/9, 0]];
    
    A12 = [[0, 0, 0, 0];...
           [1, 0, 0, 0];...
           [1/27, 8/27, 0, 0]];
    
    A21 = [[0, 0, 0];...
           [3/4, 0, 0];...
           [-2/9, 7/18, 0];...
           [-958967056548636808/9589670565486368079, 407090074900625/1274461328324957, 234959228487233/688077905161868]] ;
    
    A22 = [[0, 0, 0, 0]; ...
           [3/4, 0, 0, 0];...
           [-19/54, 14/27, 0, 0];...
           [29395718240647530075/293957182406475300749, 710656243935571/1227645422423387, -103025043292069/873209629914535, 0]];
    
    A = {A11, A12; A21, A22};


    % Coupling happens through the A coefficients
    G11 = [[1/3, 0, 0];...
           [-7/12, 1/2, 0]; ...
           [1/36, -1/6, 1/2]];
    
    G12 = [[0, 0, 0, 0];...
           [0, 0, 0, 0];...
           [0, 0, 0, 0]];
    
    G21 = [[0, 0, 0];...
           [0, 0, 0];...
           [0, 0, 0];...
           [0, 0, 0]];
    
    G22 = [[1/3, 0, 0, 0];...
           [-11/24, 1/2, 0, 0];...
           [5/12, -7/18, 1/2, 0];...
           [279874718980825/837824805432953, -303920103374805/670268586942818, ...
           -1002804453331247650/10028044533312476499, 156699807095614/247313771775319]];
    
    G = {G11, G12; G21, G22};


    B1 =  [0, 1/4, 3/4];
    B2 =  [0, 4/7, 3/7, 0];
    B = {B1,B2};
    
    BH1 =  [0, 1/4, 3/4];
    BH2 =  [65014523725041/1096981956944407, 581537073614995/1116665009524167, 144930497064493/452974281378310, 7203692236327067202/72036922363270672021];
    BH  = {BH1,BH2};
    
    % Size of problem
    N = length(y0);
    
    % Stage vectors
    K_s = zeros(N, max(s(:)), ps);
    U_s = zeros(N, max(s(:)), ps);
    
    % create an array of function/matrix handles
    rhsFuns = {rhsFun1, rhsFun2};
    jacFuns = {jacFun1, jacFun2};

    % optional arguments
    % Set the minimum number of basis vectors for the method
    if(~exist('MBasisVectors','var'))
        MBasisVectors = 1;
    end
    
    U_s(:, 1, 1) = y0;   % partition 1
    U_s(:, 1, 2) = y0;   % partition 2
    
    % krylov steps
    krySteps = 0;
    
    % phi_count
    phiCount = 0;
    
    % Jacobians or their matrix-vector product functions
    % for all partitions.
    JacJ = cell(ps,1);
    
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
    for i = 1 : max(s(:))
        % Do for each partition
        % First Exponential
        % TODO: This can be refactored to accept an array of
        % function handles and use that instead
        for j = 1:ps
            
            if i <= s(j)
                K_s(:, i , j) = 0;
                for m = 1:ps
                    for l = 1:min(i-1, s(m))
                        K_s(:, i , j) = K_s(:, i , j) + G{j,m}(i, l) * K_s(:, l, m);
                    end
                end
                
                if (~MatrixFree)
                    K_s(:, i , j) = dt * JacJ{j} *  K_s(:, i , j);
                else
                    K_s(:, i , j) = dt * JacJ{j}(K_s(:, i , j));
                end
                
                K_s(:, i , j) = K_s(:, i , j) + dt * rhsFuns{j}(U_s(:, i, j));
                
                % Compute the Krylov basis matrices for column the K's
                [V, H, M] = ArnoldiAdapt(JacJ{j}, K_s(:, i , j), N, dt * G{j,j}(i,i), MatrixFree, ...
                    NBasisVectors, relTol, MBasisVectors, adaptiveKrylov, symmjac);
                
                krySteps = krySteps + M^2;
                phiCount = phiCount + 1;
                
                % Compute Phi_1
                K_s(:, i , j) = Phi_1(G{j, j}(i, i), K_s(:, i , j), dt, V, H, M);
            end
        end
        
        % The two loops (above and below) need to be
        % kept independent. The loop below is dependent
        % on all of the computations performed in the
        % loop above.
        for j = 1:ps
            if i + 1 <= s(j)
                U_s(:, i + 1, j) =  y0;
                for m = 1:ps
                    for l = 1:min(i, s(m))
                        U_s(:, i + 1, j) = U_s(:, i + 1, j) + A{j, m}(i + 1, l) * K_s(:, l, m);
                    end
                end
            end
        end
    end
    
    y    = y0;
    yerr = zeros(size(y0));
    for m = 1:ps
        for l = 1:s(m)
            y    = y + B{m}(l) * K_s(:, l, m);
            yerr = yerr + (B{m}(l) - BH{m}(l)) * K_s(:, l, m);
        end
    end
    
    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/phiCount;
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
    
    % if v is zer0, return 0
    if(norm(v) == 0)
        phi1 = v;
        return
    end
    
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
