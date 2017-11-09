% 5th order Runge-Kutta exponential propagation integrator.
% Based on:
% [1] M. Tokman, Efficient integration of large stiff systems of ODEs with
%     exponential propagation iterative (EPI) methods.
%     J. Comput. Phys. 213 (2006) 748-776
%
% Parameters:
%   f         function handle to the rhs function
%   j         function handle to the Jacobian function
%   y0        initial conditions
%   tspan     2-element vector for the beginning and end time
%   h    step size
%   maxBasis  maximum basis size allowed for the Arnoldi iteration
%   varargin  optional parameters passed to f and j functions
%
% Returns:
%   y         final data at final time as evaluated by the 4th order method
%
% Evaluates phi-functions using Highahm high-order Pade approximation of
% the exponential function by using augmented H matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above parameters are not applicable to the parameters of the method below
% changes were made to the original code to work with the MATLODE framework.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, yerr, ISTATUS] = epirk5P1bVarSingleStep(y0, dt, rhsFun, jacFun, f, ...
                                                     MatrixFree, NBasisVectors, ISTATUS, absTol, relTol, adaptiveKrylov, symmjac, MBasisVectors, varargin)
    % Coefficients.
    a11 = 0.3512959269505819309205666343669703185949150765418610822442;
    a21 = 0.8440547201165712629778740943954345394870608397760182695023;
    a22 = 1.6905891609568963623795529686906583866077099749695686415193;

    b1 = 1.0000000000000000000000000000413981626492407621654383185164;
    b2 = 1.2727127317356892397444638436834144153459435577242007287745;
    b3 = 2.2714599265422622275490331068926657550105500849448098122502;

    g11 = 0.3512959269505819309205666513780525924758396470877807385413;
    g21 = 0.8440547201165712629778740927319685574987844036261402099979;
    g22 = 1.0;
    g31 = 1.0000000000000000000;
    g32 = 0.7111109536436687035914924546923329396898123318411132877913;
    g33 = 0.6237811195337149480885036224182618420517566408959269219908;
    g32b = 0.5;
    g33b = 1.0;
    %g32b = 0.640284661805129;
    %g33b = 0.75;
    %g32b = 0.584170797083078;
    %g33b = 0.85;

    % Main integration loop.
    N = length(y0);

    zeroVec = zeros(N,1);

    h = dt;

    hF = h * f;  % calculate right hand side - MATLODE already computes it.

    if(~MatrixFree)
      % Method wants jacFun to evaluate to the full jacobian, not calculate jacobian vector products however.
      hA = h * jacFun;  % jacobian is already computed and passed in
    else
      error('epirk5P1bVarSingleStep needs full Jacobian and is not matrix free.');
    end

    kryTol = 0.1 * min(absTol, relTol * norm(y0)) / h;
    %kryTol = min(absTol, relTol * norm(y0));
    krySteps = 0;
    
    % Stage 1.
    [temp1, stats] = phipmPaul([g11 g21 g31] , hA, [zeroVec, hF], kryTol, false, 1);
    r1 = y0 + a11*temp1(:,1);
    %    krySteps = max(stats(3), krySteps);
    krySteps = krySteps + stats(3)^2;

    % Stage 2.
    hb1 = h * feval(rhsFun, r1, varargin{:}) - hF - hA*(r1 - y0);
    [temp2, stats] = phipmPaul([g32b g32 g22] ,hA, [zeroVec, hb1], kryTol, false, 1);
    r2 = y0 + a21*temp1(:,2) + a22*temp2(:,3);
    %    krySteps = max(stats(3), krySteps);
    krySteps = krySteps + stats(3)^2;

    % Calculate next time step.
    hb2 = -2*hb1 + h * feval(rhsFun, r2, varargin{:}) - hF - hA*(r2 - y0);
    [temp3, stats] = phipmPaul([g33 g33b], hA, [zeroVec, zeroVec, zeroVec, hb2], kryTol, false, 1);
    fifthOrder = b2*temp2(:,2) + b3*temp3(:,1);
    y = y0 + b1*temp1(:,3) + fifthOrder;
    %    krySteps = max(stats(3), krySteps);
    krySteps = krySteps + stats(3)^2;

    ISTATUS.Nkdim = ISTATUS.Nkdim + krySteps/3; % Add number of
                                                % krylov steps

    % Estimate error.
    yerr = fifthOrder - b2*temp2(:,1) - b3*temp3(:,2);
end

function z = Psi1(A)
    N = length(A);
    e1 = zeros(N,1);
    e1(1,1) = 1;
    A(1:N,N+1) = e1;
    A(N+1,1:N+1) = zeros(1,N+1);
    temp = expmhigham(A);
    z = temp(1:N,end);
end

function z = Psi2(A)
    N = length(A);
    e1 = zeros(N,1);
    e1(1,1) = 1;
    A(1:N,N+1) = e1;
    A(N+1,1:N+1) = zeros(1,N+1);
    temp = expmhigham(A);
    z = temp(1:N,end);
end

function z = Psi3(A)
    N = length(A);
    e1 = zeros(N,1);
    e1(1,1) = 1;
    A(1:N,N+1) = e1;
    A(N+1:N+2,N+2:N+3) = eye(2);
    A(N+3,1:N+3) = zeros(1,N+3);
    temp = expmhigham(A);
    z = temp(1:N,end);
end

% Helper function to calculate the matrix functions, e.g. phi1, phi2, etc.
% Assumes final two terms are phin(H) * e1, n = 0,1,2.
function v = CalcMatrixFunc(f, V, H, b)
    v = norm(b) * V * feval(f, H);
end
