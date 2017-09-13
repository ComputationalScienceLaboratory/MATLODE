%% erow4AdaptMatFree
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
function [y, yerr, ISTATUS] = erow4SingleStep(y0, dt, rhsFun, J, ...
                                              f, MatrixFree, NBasisVectors, ISTATUS, absTol, relTol)
    % Based on the form for EROW4 given in "COMPARATIVE PERFORMANCE OF
    % EXPONENTIAL, IMPLICIT, AND EXPLICIT INTEGRATORS FOR STIFF SYSTEMS
    % OF ODES"  by Loffeld and Tokman.  Derived in "Exponential Rosenbrock-
    % type methods" by Hochbruck Ostermann and Schweitzer.
    N = length(y0);
    y = y0;

    normf = norm(f,2);                          % norm of f(y0)

    %----------------------------------------------------------------------
    % Set poor Tolerance level for Arnoldi
    arnoldiTol = 1e-12;
    %----------------------------------------------------------------------

    krySteps = 0;
    
    % Compute the Krylov basis matrices
    [V0, H0, M0] = ArnoldiAdapt(J, f, N, dt, MatrixFree, NBasisVectors, arnoldiTol);
    e0 = eye(M0);
    e01 = e0(:,1);
    krySteps = krySteps + M0^2;


    Hbar = [1/2*dt*H0 e01; zeros(1, M0+1)];
    tempExp = expm(Hbar);
    phi1f = normf*V0*tempExp(1:M0,end);

    Y1 = y + 1/2*dt*phi1f;

    R1 = Residual(f, rhsFun(Y1), y, Y1, J, MatrixFree);
    normR1 = norm(R1,2);
    
    [V1, H1, M1] = ArnoldiAdapt(J, R1, N, dt, MatrixFree, NBasisVectors, arnoldiTol);
    e1 = eye(M1);
    e11 = e1(:,1);
    krySteps = krySteps + M1^2;

    Hbar = [dt*H1 e11; zeros(1, M1+1)];
    tempExp = expm(Hbar);
    phi1R1 = normR1*V1*tempExp(1:M1,end);

    Hbar = [dt*H0 e01; zeros(1, M0+1)];
    tempExp = expm(Hbar);
    phi1f = normf*V0*tempExp(1:M0,end);

    Y2 = y + dt*phi1f + dt*phi1R1;

    R2 = Residual(f, rhsFun(Y2), y, Y2, J, MatrixFree);
    normR2 = norm(R2,2);

    [V2, H2, M2] = ArnoldiAdapt(J, R2, N, dt, MatrixFree, NBasisVectors, arnoldiTol);
    e2 = eye(M2);
    e21 = e2(:,1);
    krySteps = krySteps + M2^2;

    Hbar = [dt*H1 e11 zeros(M1,3); zeros(4,M1+1), [eye(3); zeros(1,3)]];
    tempExp = expm(Hbar);
    phi3R1 = normR1*V1*tempExp(1:M1, end-1);
    phi4R1 = normR1*V1*tempExp(1:M1, end);

    Hbar = [dt*H2 e21 zeros(M2,3); zeros(4,M2+1),[eye(3); zeros(1,3)]];
    tempExp = expm(Hbar);
    phi3R2 = normR2*V2*tempExp(1:M2, end-1);
    phi4R2 = normR2*V2*tempExp(1:M2, end);

    y = y + dt*phi1f + dt*(16*phi3R1 - 48*phi4R1 - 2*phi3R2 + 12*phi4R2);
    % y_hat = y + dt * (phi1f + 16 * phi3R1 - 2 * phi3R2)
    yerr = dt * (-48 * phi4R1 +12 * phi4R2);

    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/3;
end



function r = Residual(fn, fi, yn, yi, Jn, MatrixFree)
    if( ~MatrixFree )
        r = fi - fn - Jn*(yi - yn);
    else
        r = fi - fn - Jn(yi - yn);
    end
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
