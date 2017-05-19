function [ jac ] = swe_Jacobian_mex( t,y )
    
    % ATTENTION: (NX+2)*(NY+2)*3 = 3468
    % Sparse parameters will change with model parameters. This needs to be
    % generalized in the future.

     [row, col, val ] = swe_jac_mex(t,y);
     jac = sparse(row,col,val,3468,3468);

return;

