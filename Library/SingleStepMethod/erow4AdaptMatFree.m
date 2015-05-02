function [y Y1 Y2] = erow4Adapt(y0, dt, Tspan, rhsFun, jacFun)
% Based on the form for EROW4 given in "COMPARATIVE PERFORMANCE OF 
% EXPONENTIAL, IMPLICIT, AND EXPLICIT INTEGRATORS FOR STIFF SYSTEMS 
% OF ODES"  by Loffeld and Tokman.  Derived in "Exponential Rosenbrock-
% type methods" by Hochbruck Ostermann and Schweitzer.
t = [Tspan(1):dt:Tspan(end)];
N = length(y0);
y = y0;
q = Tspan(1);

for r = 1:length(t)-1
   
    f = rhsFun(0,y);
    J = @(w)jacFun(0,y,w); 
    normf = norm(f,2);                          % norm of f(y0)

    % Compute the Krylov basis matrices
    [V0 H0 M0] = ArnoldiAdaptMatFree(J, f, N, dt);
    e0 = eye(M0);
    e01 = e0(:,1);
    
    Hbar = [1/2*dt*H0 e01; zeros(1, M0+1)];
    tempExp = expm(Hbar);
    phi1f = normf*V0*tempExp(1:M0,end);
    
    Y1 = y + 1/2*dt*phi1f;
    
    R1 = Residual(f, rhsFun(0,Y1), y, Y1, J);
    normR1 = norm(R1,2);
    
    [V1 H1 M1] = ArnoldiAdaptMatFree(J, R1, N, dt);
    e1 = eye(M1);
    e11 = e1(:,1);
    
    
    Hbar = [dt*H1 e11; zeros(1, M1+1)];
    tempExp = expm(Hbar);
    phi1R1 = normR1*V1*tempExp(1:M1,end);
    
    Hbar = [dt*H0 e01; zeros(1, M0+1)];
    tempExp = expm(Hbar);
    phi1f = normf*V0*tempExp(1:M0,end);
    
    Y2 = y + dt*phi1f + dt*phi1R1;
    
    R2 = Residual(f, rhsFun(0,Y2), y, Y2, J);
    normR2 = norm(R2,2);
    
    [V2 H2 M2] = ArnoldiAdaptMatFree(J, R2, N, dt);
    e2 = eye(M2);
    e21 = e2(:,1);
  
    Hbar = [dt*H1 e11 zeros(M1,3); zeros(4,M1+1), [eye(3); zeros(1,3)]];
    tempExp = expm(Hbar);
    phi3R1 = normR1*V1*tempExp(1:M1, end-1);
    phi4R1 = normR1*V1*tempExp(1:M1, end);
    
    Hbar = [dt*H2 e21 zeros(M2,3); zeros(4,M2+1),[eye(3); zeros(1,3)]];
    tempExp = expm(Hbar);
    phi3R2 = normR2*V2*tempExp(1:M2, end-1);
    phi4R2 = normR2*V2*tempExp(1:M2, end);
    
    y = y + dt*phi1f + dt*(16*phi3R1 - 48*phi4R1 - 2*phi3R2 + 12*phi4R2);
    
end

return



function r = Residual(fn, fi, yn, yi, Jn);

r = fi - fn - Jn(yi - yn);

return