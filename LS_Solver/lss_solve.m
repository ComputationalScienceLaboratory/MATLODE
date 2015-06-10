%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ ISING X ] = lss_solve( trans, rhs, e, ISING )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Solving the system (HGamma-Jac)*X=RHS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   X = e\rhs;
   %Xout = X';
   
return;

