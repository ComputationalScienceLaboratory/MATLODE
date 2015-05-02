function [y yerr] = exp4SingleStep(y, dt, rhsFun, J, ...
                            f, MatrixFree, NBasisVectors)

% Implements exp4 given in equation (5.8) Hochbruck:  "Exponential 
% Integrators For Large Systems of Differential Equations".


N = length(y);

k = zeros(N, 7);
d = zeros(N,2);
w = zeros(N,2);
u = zeros(N,2);

normf = norm(f,2);                          % norm of f(y0)
% Compute the Krylov basis matrices
[V, H, M] = ArnoldiAdapt(J, f, N, dt, MatrixFree, NBasisVectors);
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

if( ~MatrixFree )
    d(:,1) = rhsFun(u(:,1)) - f - dt*J*w(:,1);
else
    d(:,1) = rhsFun(u(:,1)) - f - dt*J(w(:,1));
end

normd = norm(d(:,1));

%  Compute new Krylov basis matrices
[V, H, M] = ArnoldiAdapt(J, d(:,1), N, dt, MatrixFree, NBasisVectors);

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

if( ~MatrixFree )
    d(:,2) = rhsFun(u(:,2)) - f - dt*J*w(:,2);
else
    d(:,2) = rhsFun(u(:,2)) - f - dt*J(w(:,2));
end

normd = norm(d(:,2));
%  Compute new Krylov basis matrices
[V, H, M] = ArnoldiAdapt(J, d(:,2), N, 1/3*dt, MatrixFree, NBasisVectors);

e1 = [1; zeros(M-1, 1)];                    % first canonical basis vector
Hbar = [1/3*dt*H e1; zeros(1, M+1)];
tempExp = expm(Hbar);
k(:,7) = normd*V*tempExp(1:M, end);

y = y + dt*( k(:,3) + k(:,4) - 4/3*k(:,5) + k(:,6) + 1/6*k(:,7) );
yerr = y + dt*( k(:,3) - 1/2*k(:,4) - 2/3*k(:,5) + 1/2*k(:,6) + 1/2*k(:,7));
yerr = y - yerr;





