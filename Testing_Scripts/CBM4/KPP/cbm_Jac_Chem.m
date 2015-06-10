
% Wrapper for calling the sparse ODE Jacobian routine
% in a format required by Matlab's ODE integrators

function J =  cbm_Jac_Chem(T, Y)   
  
  global TIME FIX RCONST 
%  To call the mex file uncomment one of the following lines:
%     1) LU prefix if SPARSE_LU_ROW option was used in code generation
%  global LU_IROW LU_ICOL 
%     2) if SPARSE_ROW option was used in code generation
%  global IROW ICOL 
  
  Told = TIME;
  TIME = T;
  cbm_Update_SUN;
  cbm_Update_RCONST;
  
%  This line calls the Matlab ODE Jacobian routine  
  J = cbm_Jac_SP( Y, FIX, RCONST );
  
%  To call the mex routine instead, comment the line above and uncomment one of the following lines:
%     1) LU prefix if SPARSE_LU_ROW option was used in code generation
%  J = sparse( LU_IROW, LU_ICOL, ...
%        cbm_mex_Jac_SP( Y, FIX, RCONST ), 32, 32); 
%     2) if SPARSE_ROW option was used in code generation
%  J = sparse( IROW, ICOL, ...
%        cbm_mex_Jac_SP( Y, FIX, RCONST ), 32, 32); 

  TIME = Told;
  
return              
