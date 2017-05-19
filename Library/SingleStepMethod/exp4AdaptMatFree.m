%% exp4AdaptMatFree
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
function [y ksize] = exp4AdaptMatFree(y0, dt, Tspan, rhsFun, jacFun)

% Implements exp4 given in equation (5.8) Hochbruck:  "Exponential 
% Integrators For Large Systems of Differential Equations".


t = [Tspan(1):dt:Tspan(end)];
N = length(y0);
y = y0;
q = Tspan(1);

for r = 1:length(t)-1
    
    k = zeros(N, 7);
    d = zeros(N,2);
    w = zeros(N,2);
    u = zeros(N,2);

    f = rhsFun(0,y);                           % RHS evaluation
%     J = jacFun(0,y);                           % Jacobian evaluation
    J = @(w)jacFun(0,y,w);
    dFdT = zeros(N, 1);                         % problem is autonomous
    normf = norm(f,2);                          % norm of f(y0)

    % Compute the Krylov basis matrices
    [V H M] = ArnoldiAdaptMatFree(J, f, N, dt);
    ksize(1,r) = M;
    % Constructed using an augmented H as in Hochbruck "Exponential 
    % Integrators For Large Systems of Differential Equations" page 17.
    e1 = [1; zeros(M-1, 1)];                    % first canonical basis vector
    Hbar = [1/3*dt*H e1; zeros(1, M+1)];
    tempExp = expm(Hbar);
    k(:,1) = normf*V*tempExp(1:M,end);

    Hbar = [2/3*dt*H e1; zeros(1,M+1)];
    tempExp = expm(Hbar);
    k(:,2) = normf*V*tempExp(1:M, end);

    Hbar = [dt*H e1; zeros(1,M+1)];
    tempExp = expm(Hbar);
    k(:,3) = normf*V*tempExp(1:M, end);

    w(:,1) = -7/300*k(:,1) + 97/150*k(:,2) - 37/300*k(:,3);

    u(:,1) = y + dt*w(:,1);

    d(:,1) = rhsFun(0, u(:,1)) - f - dt*J(w(:,1));
    
    normd = norm(d(:,1));

    %  Compute new Krylov basis matrices
    [V H M] = ArnoldiAdaptMatFree(J, d(:,1), N, dt);
    ksize(2,r) = M;
    
    e1 = [1; zeros(M-1, 1)];                    % first canonical basis vector
    Hbar = [1/3*dt*H e1; zeros(1, M+1)];
    tempExp = expm(Hbar);
    k(:,4) = normd*V*tempExp(1:M, end);

    Hbar = [2/3*dt*H e1; zeros(1,M+1)];
    tempExp = expm(Hbar);
    k(:,5) = normd*V*tempExp(1:M, end);

    Hbar = [dt*H e1; zeros(1,M+1)];
    tempExp = expm(Hbar);
    k(:,6) = normd*V*tempExp(1:M, end);

    w(:,2) = 59/300*k(:,1) - 7/75*k(:,2) + 269/300*k(:,3) + 2/3*( ...
                                                k(:,4) + k(:,5) + k(:,6) );

    u(:,2) = y + dt*w(:,2);

    d(:,2) = rhsFun(0,u(:,2)) - f - dt*J(w(:,2));
    
    normd = norm(d(:,2));
    %  Compute new Krylov basis matrices
    [V H M] = ArnoldiAdaptMatFree(J, d(:,2), N, 1/3*dt);
    ksize(3,r) = M;
    
    e1 = [1; zeros(M-1, 1)];                    % first canonical basis vector
    Hbar = [1/3*dt*H e1; zeros(1, M+1)];
    tempExp = expm(Hbar);
    k(:,7) = normd*V*tempExp(1:M, end);
    
    y = y + dt*( k(:,3) + k(:,4) - 4/3*k(:,5) + k(:,6) + 1/6*k(:,7) );
end

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
