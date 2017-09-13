function [psiprod, mm] = evaluatepsiprod(n, m, i, j, v, dt, g, p, jac, tol)
    c = max(g(:,j)) * dt;
    [V, H, mm] = ArnoldiAdapt(jac, v, n, c, 0, 0, tol);
    psiprod = Psi(i, j, v, n, dt, g, p, V, H, mm);
end 