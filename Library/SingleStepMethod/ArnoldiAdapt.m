%% ArnoldiAdapt
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
function [V, H, i] = ArnoldiAdapt(J, f, N, c, MatrixFree, NBasisVectors, Tol, MinBasisVectors)
%
% NBasisVectors = Max number of basis vectors to be used 
% MinBasisVectors = Min number of basis vectors to be used 
% If MinBaisVectors is not specified and NBasisVectors is 0, then
% its only dependent on the residual.
%

if(~exist('Tol','var'))
    tol = 1e-12;
else
    tol = Tol;
end

if(~exist('MinBasisVectors','var'))
    MBasisVectors = 1;
else
    MBasisVectors = MinBasisVectors;
end

testIndex = [1,2,3,4,6,8,11,15,20,27,36,46,57,70,85,100];

% K-methods with max basis size
if(NBasisVectors ~= 0)
    % Preallocate H of that size plus one as is used in lines 64,65
    H = zeros(NBasisVectors+1);
    % Preallocate V of that size plus one as is used in lines 64,65
    V = zeros(N, NBasisVectors+1);
    % Store H_Size
    H_size = NBasisVectors;
else
    % Preallocate H of a reasonably big size and increase size to next
    % testIndex size if Tol check fails. Use a size from testIndex plus
    % one as is used in lines 64,65
    H_size = 15;
    H = zeros(H_size + 1);
    V = zeros(N, H_size + 1);
end


beta = norm(f);
V(:, 1) = f/beta;

%     Arnoldi iteration
for i = 1:N
    if( ~MatrixFree )
        w = J*V(:,i);
    else
        w = J(V(:,i));
    end
    
    for j = 1:i
        H(j,i) = w'*V(:,j);
        w = w - H(j,i)*V(:,j);
    end
    H(i+1,i) = norm(w);
    V(:,i+1) = w/H(i+1,i);
    
    % Perform tolerance check only when the number of basis is not
    % set and only for the chosen testIndex values
    % With line 89, 90 the i>100 condition is redundant
    if  (NBasisVectors == 0) && (~isempty(find(testIndex == i, 1)))
        % First compute the residual
        e1 = [1; zeros(i-1, 1)]; 
        Hbar = [c*H(1:i,1:i) e1; zeros(1, i+1)];
        tempExp = expm(Hbar);
        residual = c*beta*H(i+1,i)*tempExp(i, end)/N;
        % Check if the residual is below prescribed value of tolerance
        if (residual < tol) 
            break
        else
            % Check that we are not increasing the basis size beyond NBasisVectors
            if(i == NBasisVectors)
                break;
            end
            
            % Now check if we have already come up previously allocated size
            % If so increase the size
            if(i == H_size)
                % First verify whether it is in the list of testIndex
                iFindTestIndex = find(testIndex == i);
                if(length(iFindTestIndex) == 1 && iFindTestIndex < 16)
                    % If it is in testIndex, then increase it to next
                    % bigger size if one is available
                    H_size = testIndex(iFindTestIndex + 1); 
                else
                    % Increase it by a fixed percentage: Say 10%
                    H_size = int32((H_size * 110)/100);
                    % Also add it dynamically to testIndex to ensure that
                    % residual check only happens at that index
                    testIndex(end + 1) = H_size;
                end
                % H_size plus one as is used in lines 64,65
                H_new = zeros(H_size + 1);
                H_new(1:i+1, 1:i) = H(1:i+1,1:i);
                H = H_new;
                V_new = zeros(N, H_size + 1);
                V_new(:, 1:i+1) = V(:,1:i+1);
                V = V_new;
            end
        end
    end

    if i == NBasisVectors || i == 1000
        break;
    end
end
H = H(1:i,1:i);
V = V(1:N,1:i);

return

%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt.edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 
%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>
