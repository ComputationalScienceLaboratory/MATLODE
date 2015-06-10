%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ ISING, e1 ] = lss_decomp( NVAR, hgamma, fjac  )
% Based off of lapack_decomp( hgamma, ising ) [LS.Solver.F90]
        
    e1 = -fjac + eye(NVAR,NVAR)*hgamma;
        
    condNum = rcond(e1);
    if ( condNum <= eps )
        ISING = true; % ising -> 1
    else
        ISING = false; % ising -> 0
    end
   

return;

