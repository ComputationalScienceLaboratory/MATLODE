% Originally highPhi by Paul
% This is method uses the recursive definition of 
% Phi functions to build up various Phi_Ks for 
% different K values.
% The recursive definition is given in 
% i) Exponential-Krylov methods for ODEs, P. Tranquilli and A. Sandu (under
% Lemma 1)
% ii) Exponential Rosenbrock-Type Methods, M. Hochbruck, A. Ostermann
% and J. Schweitzer (Equation 2.5)
% Note: This method will sort the phiNums and return the sorted phiNums 
% along with the phiKs.
function [phiNums, phi_ks] = phiKs(A, M, phiNums, LA, UA, pA, dtA, gA)
    index = 1;
    phiNums = sort(phiNums);
    max_k = phiNums(end);

    % Pre-allocate memory.
    phi_ks = zeros(M, M,length(phiNums));

    % Assuming phiNums >= 1
    for i = 1:max_k
        % Compute the phi function for k = 1
        if i == 1
            if ~exist('LA', 'var') || ~exist('UA', 'var') || ...
                    ~exist('pA', 'var')
                phiNew = phiFun(A);
            else
                phiNew = phiFun(A, LA, UA, pA, dtA, gA);
            end
        else
            if ~exist('LA', 'var') || ~exist('UA', 'var') ...
                || ~exist('pA', 'var') || ~exist('dtA','var') || ~exist('gA','var')
                phiOld0 = 1/factorial(i-1)*eye(M);
                phiNew = (phiOld - phiOld0)/A;
            else
                phiOld0 = 1/factorial(i-1)*eye(M);
                phiDiff = (phiOld - phiOld0);
                phiNew = UA\(LA\(phiDiff(pA, :)/(dtA*gA)));            
            end
        end
        
        if  i == phiNums(index) 
            phi_ks(:,:,index) = phiNew;
            index = index + 1;
        end
        
        phiOld = phiNew;
    end
return