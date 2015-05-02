function [ ISING X ] = lss_solve( trans, rhs, e, ISING )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Solving the system (HGamma-Jac)*X=RHS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   X = e\rhs;
   %Xout = X';
   
return;

