function [L, U, p, sing] = lu_decomposition(A)

    N = size(A,1);
    [L, U, p] = lu(A, 'vector');
    sing = (nnz(abs(diag(L)) > 2*eps) ~= N || nnz(abs(diag(U)) > 2*eps) ~= N);

return;

