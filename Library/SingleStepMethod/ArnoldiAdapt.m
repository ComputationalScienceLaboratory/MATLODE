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
function [V H i] = ArnoldiAdapt(J, f, N, c, MatrixFree, NBasisVectors)


beta = norm(f);
V(:, 1) = f/beta;
tol = 1e-12;
testIndex = [1,2,3,4,6,8,11,15,20,27,36,46,57,70,85,100];
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
    
    if min(abs(testIndex - j)) == 0 || (j > 100) 
        e1 = [1; zeros(i-1, 1)]; 
        Hbar = [c*H(1:i,1:i) e1; zeros(1, i+1)];
        tempExp = expm(Hbar);
        residual = c*beta*H(i+1,i)*tempExp(i, end)/N;
        if residual < tol || j == NBasisVectors
            break
        end
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
